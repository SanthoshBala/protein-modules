package org.hupo.psi.mi.psicquic.view.webapp.cytoscape;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.view.webapp.util.UrlUtils;
import psidev.psi.mi.tab.PsimiTabException;
import psidev.psi.mi.tab.converter.tab2graphml.Tab2Cytoscapeweb;

import javax.servlet.ServletException;
import javax.servlet.ServletOutputStream;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.HttpServlet;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * Servlet enabling the conversion of a MITAB Stream into GraphML.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 1.2
 */
public class GraphmlServlet extends HttpServlet {

	private static final Log log = LogFactory.getLog(GraphmlServlet.class);

	@Override
	public void service(ServletRequest request, ServletResponse response) throws ServletException, IOException {

		String mitabUrl = request.getParameter("mitabUrl");

		String url = mitabUrl.substring(0, mitabUrl.lastIndexOf("/"));
		String query = mitabUrl.substring(mitabUrl.lastIndexOf("/"));

		if (query != null) {
			query = UrlUtils.encodeUrl(query);

			String encodedUrl = url + query;

			if (log.isDebugEnabled()) log.debug("encodedUrl: " + encodedUrl);

			final URL inputUrl = new URL(encodedUrl);
			final InputStream is = inputUrl.openStream();

			ServletOutputStream stream = response.getOutputStream();
			try {
				response.setContentType("text/plain");

				final Tab2Cytoscapeweb tab2Cytoscapeweb = new Tab2Cytoscapeweb();

				String output = null;
				try {
					output = tab2Cytoscapeweb.convert(is);
				} catch (PsimiTabException e) {
					throw new IllegalStateException("Could not parse input MITAB.", e);
				}
				stream.write(output.getBytes());
				stream.flush();
			} finally {
				is.close();
			}
		}
	}//We don't have url
}
