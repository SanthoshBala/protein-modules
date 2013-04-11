package org.hupo.psi.mi.psicquic.view.webapp.controller.servlet;

import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 03/10/2012
 * Time: 14:21
 * To change this template use File | Settings | File Templates.
 */
public class AppMainHeaderServlet extends HttpServlet {

	/**
	 * Performs an HTTP GET request
	 * @param request The {@link javax.servlet.http.HttpServletRequest}
	 * @param response The {@link javax.servlet.http.HttpServletResponse}
	 */
	public void doGet (HttpServletRequest request, HttpServletResponse response) throws IOException {
		response.setContentType("text/html");

		URL url = new URL("http://wwwdev.ebi.ac.uk/frontier/template-service/templates/services/header/main");
		// make post mode connection
		URLConnection urlc = url.openConnection();
		//urlc.setDoOutput(true);
		//urlc.setAllowUserInteraction(false);
		//OutputStreamWriter wr = new OutputStreamWriter(urlc.getOutputStream());
		//wr.flush();

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
}
