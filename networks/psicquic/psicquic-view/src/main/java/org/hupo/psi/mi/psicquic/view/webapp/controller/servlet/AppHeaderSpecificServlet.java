package org.hupo.psi.mi.psicquic.view.webapp.controller.servlet;

import javax.servlet.ServletConfig;
import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.net.URL;
import java.net.URLConnection;

public class AppHeaderSpecificServlet extends HttpServlet {

	/**
	 * Initialize the <code>AppHeaderSpecificServlet</code>
	 * @param servletConfig The Servlet configuration passed in by the servlet container
	 */
	public void init(ServletConfig servletConfig) throws ServletException {
		super.init(servletConfig);
	}

	/**
	 * Performs an HTTP GET request
	 * @param request The {@link javax.servlet.http.HttpServletRequest}
	 * @param response The {@link javax.servlet.http.HttpServletResponse}
	 */
	public void doGet (HttpServletRequest request, HttpServletResponse response) throws IOException {
		response.setContentType("text/html");

		URL url = new URL("http://wwwdev.ebi.ac.uk/frontier/template-service/templates/header/specific");
		// make post mode connection
		URLConnection urlc = url.openConnection();
		urlc.setDoOutput(true);
		urlc.setAllowUserInteraction(false);
		OutputStreamWriter wr = new OutputStreamWriter(urlc.getOutputStream());
		wr.write(getWebConfigurationJSON());
		wr.flush();

		// retrieve result
		BufferedReader br1 = new BufferedReader(new InputStreamReader(urlc.getInputStream(), "UTF-8"));
		StringBuilder sb = new StringBuilder();
		String line;
		while ((line = br1.readLine()) != null) {
			sb.append(line);
			sb.append("\n");
		}
		br1.close();
		response.getWriter().write(sb.toString());
	}

	private String getWebConfigurationJSON() throws IOException {
		StringBuilder sb = new StringBuilder();
		ServletContext servletContext = this.getServletContext();
		String pathContext = servletContext.getRealPath("WEB-INF/webconfig/config.json");
		FileInputStream fstream = new FileInputStream(pathContext);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		//Read File Line By Line
		while ((strLine = br.readLine()) != null)   {
			sb.append(strLine);
		}
		return sb.toString();
	}
}
