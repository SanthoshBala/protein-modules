package org.hupo.psi.mi.psicquic.view.webapp.io.download;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.view.webapp.io.DownloadUtils;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;

import javax.faces.application.FacesMessage;
import javax.faces.context.FacesContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.*;

public class DownloadAllServlet extends HttpServlet {


	private static final Log log = LogFactory.getLog(DownloadAllServlet.class);

	private DownloadUtils downloadAllUtils = new DownloadAllUtils();

	public static final int BUFFER_TAM = 2048;


	@Override
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		handleRequest(request, response);

	}

	@Override
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		handleRequest(request, response);

	}

	/**
	 * Process the given request, generating a response.
	 *
	 * @param request  current HTTP request
	 * @param response current HTTP response
	 * @throws javax.servlet.ServletException in case of general errors
	 * @throws java.io.IOException            in case of I/O errors
	 */
	public void handleRequest(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

//		final String serviceURL = request.getParameter("serviceURL");
		String query = request.getParameter("query");
		String format = request.getParameter("format");
		String dir = request.getParameter("dir");

//
//		if (serviceURL == null) {
//			FacesContext context = FacesContext.getCurrentInstance();
//			String msg = "The serviceURL is not provided.";
//			String details = "The file can not be downloaded.";
//
//			if (context != null) {
//				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
//				context.addMessage(null, facesMessage);
//			}
//
//			log.error(msg + " " + details);
//			return;
//		}
//
		if (query == null) {
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "The query parameter is not provided.";
			String details = "It will be used * by default.";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}

			log.warn(msg + " " + details);
			query = "*";
		}

		if (format == null) {
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "The format parameter is not provided.";
			String details = "It will be used MITAB 25 by default.";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}

			log.warn(msg + " " + details);
			format = PsicquicSimpleClient.MITAB25;
		}

		if (dir == null) {
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "The dir is not valid";
			String details = "You can not download the file";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}

			log.error(msg + " " + details);

			return;
		}

		String filename = getFileName(query);
		String extension = getExtension(format);
		String contentType = getContentType(format);

//

//		catch(Exception e){
//			FacesContext context = FacesContext.getCurrentInstance();
//			String msg = "Temporarily impossible to retrieve results";
//			String details = "Error while building results based on the query: '" + query + "'";
//
//			if (context != null) {
//				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
//				context.addMessage(null, facesMessage);
//			}
//
//			log.error(details + " SERVICE RESPONSE: " + e.getMessage());
//		}

		// count the results
		// Some JSF component library or some Filter
		// might have set some headers in the buffer beforehand.
		// We want to get rid of them, else it may collide.
		response.reset();

		if (contentType != null)
			response.setContentType(contentType);

		// double quotes are needed in case the filename is long,
		// otherwise the filename gets
		// truncated in Firefox.
		response.setHeader("Content-Disposition",
				"attachment; filename=\"" + filename + "." + extension + "\"");

		OutputStream out = null;

		try {

			out = response.getOutputStream();
			String header = writeHeaderInMitab(format);

			if (header != null) {
				out.write(header.getBytes());
				out.flush();   // Write system dependent end of line.
			}

			File folder = new File(dir);

			for (File sourceFile : folder.listFiles()) {
				if (!sourceFile.isDirectory()) {
					InputStream in = null;
					try {
						in = new FileInputStream(sourceFile);
						// Transfer bytes from in to out
						byte[] buf = new byte[BUFFER_TAM];

						while (in.read(buf) > 0) {
							out.write(buf);
						}
						out.flush();
					} finally {
						if (in != null) {
							try {
								in.close();
								if(!sourceFile.delete()){
									log.info("The file: " + sourceFile + "can not be deleted.");
								}
							} catch (IOException e) {
								//Keep it quiet
							}
						}
					}
				}
			}

			if(!folder.delete()){
				log.info("The directory: " + folder + "can not be deleted.");
			}

		} finally {
			if (out != null) {
				try {
					out.close();
				} catch (IOException e) {
					//Keep it quiet
				}
			}
		}
	}

	private String getFileName(String query) {
		return downloadAllUtils.getFileName(query);
	}

	private String getExtension(String format) {
		return downloadAllUtils.getExtension(format);
	}

	private String getContentType(String format) {
		return downloadAllUtils.getContentType(format);
	}

	private String writeHeaderInMitab(String format) {
		return downloadAllUtils.writeHeaderInMitab(format);
	}

}
