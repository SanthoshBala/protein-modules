package org.hupo.psi.mi.psicquic.view.webapp.io.download;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.myfaces.trinidad.context.RequestContext;
import org.apache.myfaces.trinidad.event.ReturnEvent;
import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.view.webapp.application.PsicquicThreadConfig;
import org.hupo.psi.mi.psicquic.view.webapp.controller.BaseController;
import org.hupo.psi.mi.psicquic.view.webapp.controller.config.DetectProxy;
import org.hupo.psi.mi.psicquic.view.webapp.controller.config.PsicquicViewConfig;
import org.hupo.psi.mi.psicquic.view.webapp.controller.search.SearchController;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;

import javax.faces.application.FacesMessage;
import javax.faces.context.FacesContext;
import javax.faces.event.ValueChangeEvent;
import javax.faces.model.SelectItem;
import java.io.*;
import java.net.Proxy;
import java.util.*;
import java.util.concurrent.*;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 30/10/2012
 * Time: 10:15
 * To change this template use File | Settings | File Templates.
 */

@Controller("downloadAllController")
public class DownloadAllController extends BaseController {

	private List<SelectItem> selectAllItems;

	private String selectedFormat;

	private static final Log log = LogFactory.getLog(DownloadAllController.class);


	@Autowired
	private SearchController searchController;

	@Autowired
	private PsicquicViewConfig config;

	private ArrayList<Future> runningTasks;


	public DownloadAllController() {
		this.selectAllItems = new ArrayList<SelectItem>();
	}

	private void createAllList() {
		ResourceBundle rb = ResourceBundle.getBundle("org.hupo.psi.mi.psicquic.Messages");

		String format = "tab25";

		this.selectAllItems.add(new SelectItem("tab25", rb.getString("tab25").trim(), rb.getString("tab25").trim()));
		this.selectAllItems.add(new SelectItem("tab26", rb.getString("tab26").trim(), rb.getString("tab26").trim(), true));
		this.selectAllItems.add(new SelectItem("tab27", rb.getString("tab27").trim(), rb.getString("tab27").trim(), true));

		setSelectedFormat(format);
	}


	public void setSelectAllItems(List<SelectItem> selectAllItems) {
		this.selectAllItems = selectAllItems;
	}

	public List<SelectItem> getSelectAllItems() {
		if (selectAllItems.isEmpty())
			createAllList();

		return selectAllItems;
	}

	public String getSelectedFormat() {
		return selectedFormat;
	}

	public void setSelectedFormat(String selectedFormat) {
		this.selectedFormat = selectedFormat;
	}

	public String allDownload() {

		String download = "all";

		if (RequestContext.getCurrentInstance() != null) {
			RequestContext.getCurrentInstance().returnFromDialog(download, null);
		}

		return null;
	}

	public void forward(ReturnEvent returnEvent) {

		if (returnEvent.getReturnValue() instanceof String) {
			String value = (String) returnEvent.getReturnValue();

			FacesContext context = FacesContext.getCurrentInstance();
			String baseURL = context.getExternalContext().getRequestContextPath();
//			String serviceURL = searchController.getServicesMap().get(searchController.getSelectedServiceName()).getRestUrl();
			String query = searchController.getUserQuery().getSearchQuery();
			String format = getSelectedFormat();

			if (value.contentEquals("all")) {

				String path = getResultsPerService(query, format);
			    if(path!=null){
					String url = baseURL + "/allDownload?" +
							"&query=" + query +
							"&format=" + format +
							"&dir=" + path;

					try {
						String encodeURL = context.getExternalContext().encodeResourceURL(url);
						context.getExternalContext().redirect(encodeURL);
					} catch (Exception e) {
						final String msg = "The content of the file can not be generated";
						final String detail = "The " + url + "for download the file is not available";
						addErrorMessage(msg, detail);
						log.error("The content of the file can not be generated: " + url, e);

					} finally {
						context.responseComplete();
					}
				}
				else{
					final String msg = "The content of the file can not be generated";
					final String detail = "The url for download the file is not available";
					addErrorMessage(msg, detail);
					log.error("The temporal dir can not be generated");
					context.responseComplete();

				}
			}
		}
	}

	private String getResultsPerService(String query, String format) {

		PsicquicThreadConfig threadConfig = (PsicquicThreadConfig) searchController.getApplicationContext().getBean("psicquicThreadConfig");

		ExecutorService executorService = threadConfig.getExecutorService();

		final File tempDir = new File(config.getDownloadAllLocationFile(), Long.toString(System.nanoTime()));
		Boolean result = tempDir.mkdir();

		if(!result){
			return null;
		}

		if (runningTasks == null) {
			runningTasks = new ArrayList<Future>();
		} else {
			runningTasks.clear();
		}

		List<ServiceType> services = searchController.getServices();
		Map<String, Boolean> serviceSelectionMap = searchController.getServiceSelectionMap();


		final Map<String, File> resultFileMap = new HashMap<String, File>();


		for (final ServiceType service : services) {
			if (service.isActive() && serviceSelectionMap.get(service.getName()) != null && serviceSelectionMap.get(service.getName())) {
				final String finalQuery = query;
				final String finalFormat = format;
				Runnable runnable = new Runnable() {
					public void run() {
						File tempFile = downloadTempFile(service, finalQuery, service.getName(), finalFormat, tempDir);
						if (tempFile != null) {
							resultFileMap.put(service.getName(), tempFile);
						}
					}
				};
				runningTasks.add(executorService.submit(runnable));
			}
		}

//		ExecutorService pool = Executors.newFixedThreadPool(29);
//
//		List<ServiceType> services = searchController.getServices();
//		Map<String, Boolean> serviceSelectionMap = searchController.getServiceSelectionMap();
//		final Map<String, File> resultFileMap = new HashMap<String, File>();
//
//
//		for (final ServiceType service : services) {
//			if (service.isActive() && serviceSelectionMap.get(service.getName()) != null && serviceSelectionMap.get(service.getName())) {
//				final String finalQuery = query;
//				final String finalFormat = format;
//				Runnable runnable = new Runnable() {
//					public void run() {
//						File tempFile = downloadTempFile(service, finalQuery, service.getName(),finalFormat);
//						resultFileMap.put(service.getName(), tempFile);
//					}
//				};
//				runningTasks.add(pool.submit(runnable));
//			}
//		}
		checkAndResumePsicquicTasks();

		return tempDir.getAbsolutePath();

	}

	private void checkAndResumePsicquicTasks() {

		for (Future f : runningTasks) {
			int threadTimeOut = 5;
			try {
				f.get(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
			} catch (InterruptedException e) {
				log.error("The psicquic task was interrupted, we cancel the task.", e);
				if (!f.isCancelled()) {
					f.cancel(true);
				}
			} catch (ExecutionException e) {
				log.error("The psicquic task could not be executed, we cancel the task.", e);
				if (!f.isCancelled()) {
					f.cancel(true);
				}
			} catch (TimeoutException e) {
				log.error("Service task stopped because of time out " + threadTimeOut + "seconds.", e);

				if (!f.isCancelled()) {
					f.cancel(true);
				}
			}
		}

		runningTasks.clear();
	}

	private static File downloadTempFile(ServiceType service, String query, String name, String format, File tempDir) {

		PsicquicSimpleClient psicquicSimpleClient;
		DetectProxy detectProxy = new DetectProxy();
		Proxy proxy = detectProxy.getProxy();

		if (proxy != null) {
			psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl(), proxy);
		} else {
			psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl());
		}
		psicquicSimpleClient.setReadTimeout(40000);

		BufferedWriter out = null;
		BufferedReader in = null;
		OutputStream outputStream = null;
		File tempFile = null;


		try {

			InputStream inputStream = psicquicSimpleClient.getByQuery(query, format);

			tempFile = File.createTempFile(name, null, tempDir);
			log.info("The temporary file" + " has been created: " + tempFile);

			outputStream = new FileOutputStream(tempFile);

			out = new BufferedWriter(new OutputStreamWriter(outputStream));
			in = new BufferedReader(new InputStreamReader(inputStream));

			String userInput;

			//TODO Read by chunks instead of lines
			while ((userInput = in.readLine()) != null) {
				out.write(userInput);
				out.newLine();   // Write system dependent end of line.
			}
			out.flush();

		} catch (IOException e) {
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "Temporarily impossible to retrieve results";
			String details = "Error while building results based on the query: '" + query + "'";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}

			//TODO Check if it is a 2.7 query and explain
			log.error("Error while building results based on the query: '"
					+ query + "'" + " SERVICE (" + service.getName() + ") RESPONSE: " + e.getMessage());
//			throw new IOException("Error while building results based on the query: '" + query + "'",e);

		} finally {
			try {
				if (in != null) {
					in.close();
				}
				if (out != null) {
					out.close();
				}
			} catch (IOException e) {
				//Keep it quiet
			}
		}

		return tempFile;
	}

	public String cancelDialog() {

		if (RequestContext.getCurrentInstance() != null) {
			RequestContext.getCurrentInstance().returnFromDialog(null, null);
		}

		return null;
	}

	public void selectedFormatChanged(ValueChangeEvent valueChangeEvent) {
		if (valueChangeEvent.getNewValue() instanceof String) {
			selectedFormat = ((String) valueChangeEvent.getNewValue());

		}
	}

}
