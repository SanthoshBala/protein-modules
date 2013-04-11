package org.hupo.psi.mi.psicquic.view.webapp.controller.clustering;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.QueryResponse;
import org.hupo.psi.mi.psicquic.clustering.ClusteringServiceException;
import org.hupo.psi.mi.psicquic.clustering.InteractionClusteringService;
import org.hupo.psi.mi.psicquic.model.PsicquicSolrServer;
import org.hupo.psi.mi.psicquic.view.webapp.io.DownloadUtils;
import org.springframework.context.ApplicationContext;
import org.springframework.stereotype.Controller;
import org.springframework.web.context.support.WebApplicationContextUtils;

import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Servlet to download the data stored in the local index holding clustered data.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 1.2
 */
@Controller
public class ClusteringDownloadServlet extends HttpServlet {

    private static final Log log = LogFactory.getLog(ClusteringDownloadServlet.class);

    private InteractionClusteringService clusteringService;
	private final DownloadUtils clusteringDownloadUtils = new ClusteringDownloadUtils();


	@Override
    public void init(ServletConfig config) throws ServletException {
        super.init(config);
        ApplicationContext ctx = WebApplicationContextUtils.getRequiredWebApplicationContext(config.getServletContext());
        clusteringService = (InteractionClusteringService) ctx.getBean( "defaultInteractionClusteringService" );
    }

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

    protected void processRequest(HttpServletRequest request, HttpServletResponse response) throws IOException {

        final String jobId = request.getParameter( "jobId" );
		if( jobId == null) {
			log.error("The jobId is not provided, the file can not be downloaded.");
			return;
		}

		String query = request.getParameter( "query" );
        if( query == null) {
			log.warn("The query parameter is not provided, it will be used * by default");
			query = "*";
        }

		String format = request.getParameter( "format" );

        if( format == null) {
			log.warn("The format parameter is not provided, it will be used MITAB 25 by default");
			format = PsicquicSolrServer.RETURN_TYPE_MITAB25;
        }

		String filename = getFileName(query);
		String extension = getExtension(format);
		String contentType = getContentType(format);

		response.reset();

		if (contentType != null){
			response.setContentType(contentType);
		}

		response.setHeader("Content-Disposition",
				"attachment; filename=\"" + filename + "." + extension + "\"");

		final OutputStream stream = response.getOutputStream();

        try {
            int current = 0;
            final int batchSize = 200;
            QueryResponse resp;

			String header = clusteringDownloadUtils.writeHeaderInMitab(format);

			if(header!=null){
				stream.write(header.getBytes());
				stream.flush();   // Write system dependent end of line.
			}

            do {
                // query local data
                resp = clusteringService.query(jobId, query, current, batchSize, format );

                // Write MITAB on the output stream
                String mitab = resp.getResultSet().getMitab();

                // the job could be found/and and some results could be found associated with this job.
                if (mitab != null){
                    stream.write(mitab.getBytes());
                    stream.flush();
                }

                current += batchSize;
            } while ( current < resp.getResultInfo().getTotalResults() );

            // convert to a list of BinaryInteraction
            if( log.isInfoEnabled() ) log.info("Total MITAB lines: " + resp.getResultInfo().getTotalResults() );

        } catch ( ClusteringServiceException e ) {
            throw new RuntimeException( "Failed to query local clustered data", e );
        } finally {
            stream.flush();
            stream.close();
        }
    }

	private String getFileName(String query) {
		return clusteringDownloadUtils.getFileName(query);
	}

	private String getExtension(String format) {
		return clusteringDownloadUtils.getExtension(format);
	}

	private String getContentType(String format) {
		return clusteringDownloadUtils.getContentType(format);
	}
}
