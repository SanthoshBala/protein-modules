package org.hupo.psi.mi.psicquic.ws.utils;

import org.hupo.psi.calimocho.tab.util.MitabDocumentDefinitionFactory;
import org.hupo.psi.calimocho.xgmml.XgmmlStreamingGrapBuilder;
import org.hupo.psi.mi.psicquic.PsicquicService;
import org.hupo.psi.mi.psicquic.RequestInfo;
import org.hupo.psi.mi.psicquic.ws.IndexBasedPsicquicService;

import javax.ws.rs.WebApplicationException;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Streaming output for XGMML
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>08/08/12</pre>
 */

public class XgmmlStreamingOutput extends PsicquicStreamingOutput {

    public XgmmlStreamingOutput(PsicquicService psicquicService, String query, int first, int maxRows) {
        super(psicquicService, query, first, maxRows);
    }

    @Override
    public void write(OutputStream outputStream) throws IOException, WebApplicationException {

        int blockSize = Math.min(IndexBasedPsicquicService.BLOCKSIZE_MAX, maxResults);
        int firstResult = this.firstResult;

        XgmmlStreamingGrapBuilder graphBuilder = null;
        try {
            graphBuilder = new XgmmlStreamingGrapBuilder("PSICQUIC", "Generated from MITAB 2.5", "http://psicquic.googlecode.com");

            graphBuilder.open(outputStream, countResults());
            RequestInfo reqInfo = new RequestInfo();
            reqInfo.setResultType("psi-mi/tab25");

            OutputStream os=outputStream;

            int totalResults = -1;
            int max = firstResult + maxResults;

            do {
                reqInfo.setFirstResult(firstResult);
                reqInfo.setBlockSize(blockSize);

                try {
                    response = psicquicService.getByQuery(query, reqInfo);
                    InputStream mitab = new ByteArrayInputStream(response.getResultSet().getMitab().getBytes());
                    try{
                        graphBuilder.writeNodesAndEdgesFromMitab(mitab, MitabDocumentDefinitionFactory.mitab25());
                    }
                    finally {
                        mitab.close();
                    }
                } catch (Exception e) {
                    throw new WebApplicationException(e);
                }

                totalResults = response.getResultInfo().getTotalResults();

                firstResult = firstResult + blockSize;

            } while (firstResult < totalResults && firstResult < max);
        } catch (JAXBException e) {
            throw new WebApplicationException(e);
        } catch (XMLStreamException e) {
            throw new WebApplicationException(e);
        }
        finally {
            if (graphBuilder != null){
                try {
                    graphBuilder.close();
                } catch (XMLStreamException e) {
                    throw new WebApplicationException(e);
                }
            }
        }
    }

}
