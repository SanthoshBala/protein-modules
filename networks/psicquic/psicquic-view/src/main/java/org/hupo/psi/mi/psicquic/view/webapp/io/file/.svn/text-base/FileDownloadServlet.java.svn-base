package org.hupo.psi.mi.psicquic.view.webapp.io.file;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.view.webapp.controller.config.DetectProxy;
import org.hupo.psi.mi.psicquic.view.webapp.io.DownloadUtils;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;

import javax.faces.application.FacesMessage;
import javax.faces.context.FacesContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.net.Proxy;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 29/08/2012
 * Time: 15:28
 * To change this template use File | Settings | File Templates.
 */
public class FileDownloadServlet extends HttpServlet {

	private static final Log log = LogFactory.getLog(FileDownloadServlet.class);


	private final DownloadUtils fileDownloadUtils = new FileDownloadUtils();

	@Override
	protected void doGet(
			HttpServletRequest aRequest, HttpServletResponse aResponse
	) throws ServletException, IOException {
		processRequest(aRequest, aResponse);
	}

	@Override
	protected void doPost(
			HttpServletRequest aRequest, HttpServletResponse aResponse
	) throws ServletException, IOException {
		processRequest(aRequest, aResponse);
	}

	private void processRequest(HttpServletRequest request, HttpServletResponse response) throws IOException {

		final String serviceURL = request.getParameter("serviceURL");
		String query = request.getParameter("query");
		String format = request.getParameter("format");

		if (serviceURL == null) {
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "The serviceURL is not provided.";
			String details = "The file can not be downloaded.";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}

			log.error(msg + " " + details);
			return;
		}

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
		String filename = getFileName(query);
		String extension = getExtension(format);
		String contentType = getContentType(format);

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

		OutputStream outputStream = null;
		BufferedWriter out = null;
		BufferedReader in = null;

		try {
			outputStream = response.getOutputStream();

			// Now you can write the InputStream of the
			// file to the above OutputStream the usual way.

			PsicquicSimpleClient psicquicSimpleClient;
			DetectProxy detectProxy = new DetectProxy();
			Proxy proxy = detectProxy.getProxy();

			if(proxy != null){
				psicquicSimpleClient = new PsicquicSimpleClient(serviceURL, proxy);
			} else {
				psicquicSimpleClient = new PsicquicSimpleClient(serviceURL);
			}

			psicquicSimpleClient.setReadTimeout(40000);
			InputStream inputStream = psicquicSimpleClient.getByQuery(query, format);

			String userInput;

			out = new BufferedWriter(new OutputStreamWriter(outputStream));
			in = new BufferedReader(new InputStreamReader(inputStream));

			String header = writeHeaderInMitab(format);

			if(header!=null){
				out.write(header);
				out.newLine();   // Write system dependent end of line.
			}

			//TODO Read by chunks instead of lines
			while ((userInput = in.readLine()) != null) {
				out.write(userInput);
				out.newLine();   // Write system dependent end of line.
			}
			out.flush();
		}
		catch(Exception e){
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "Temporarily impossible to retrieve results";
			String details = "Error while building results based on the query: '" + query + "'";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}

			log.error(details + " SERVICE RESPONSE: " + e.getMessage());
			throw new IOException("Error while building results based on the query: '" + query + "'",e);
		}
		finally {
			if (in != null) {
				in.close();
			}
			if (out != null) {
				out.close();
			}
		}

	}

	private String getFileName(String query) {
		return fileDownloadUtils.getFileName(query);
	}

	private String getExtension(String format) {
		return fileDownloadUtils.getExtension(format);
	}

	private String getContentType(String format) {
		return fileDownloadUtils.getContentType(format);
	}

	private String writeHeaderInMitab(String format) {
		return fileDownloadUtils.writeHeaderInMitab(format);
	}
}
