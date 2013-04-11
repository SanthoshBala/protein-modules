package org.hupo.psi.mi.psicquic.indexing.batch;

import org.apache.commons.io.FileUtils;
import org.hupo.psi.mi.psicquic.model.server.SolrJettyRunner;
import org.junit.After;
import org.junit.Before;

/**
 * Abstract class for testing solr server
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>16/07/12</pre>
 */

public abstract class AbstractSolrServerTest {

    protected SolrJettyRunner solrJettyRunner;

    @Before
    public void startJetty() throws Exception {
        // Start a jetty server to host the solr index
        this.solrJettyRunner = new SolrJettyRunner();
        solrJettyRunner.start();
    }

    @After
    public void stopJetty() throws Exception {
        // shutdown solrJetty
        solrJettyRunner.stop();

        FileUtils.deleteQuietly(solrJettyRunner.getSolrHome());
    }
}
