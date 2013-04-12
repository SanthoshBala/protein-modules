package org.hupo.psi.mi.psicquic.model;

import org.apache.solr.client.solrj.SolrQuery;
import org.apache.solr.client.solrj.SolrServer;
import org.apache.solr.client.solrj.SolrServerException;
import org.apache.solr.common.SolrDocumentList;
import org.hupo.psi.mi.psicquic.indexing.batch.AbstractSolrServerTest;
import org.hupo.psi.mi.psicquic.indexing.batch.SolrMitabIndexer;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.batch.core.JobParametersInvalidException;
import org.springframework.batch.core.repository.JobExecutionAlreadyRunningException;
import org.springframework.batch.core.repository.JobInstanceAlreadyCompleteException;
import org.springframework.batch.core.repository.JobRestartException;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;
import psidev.psi.mi.tab.PsimiTabException;
import psidev.psi.mi.tab.PsimiTabReader;

import java.io.IOException;
import java.util.Collection;

/**
 * Tester of MitabInputStream
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/07/12</pre>
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration( locations = {"classpath*:/META-INF/psicquic-spring.xml",
        "classpath*:/jobs/psicquic-indexing-spring-test.xml"})
public class MitabInputStreamTest  extends AbstractSolrServerTest {

    @Autowired
    private SolrMitabIndexer solrMitabIndexer;
    private PsimiTabReader mitabReader = new PsimiTabReader();

    @Test
    public void test_create_mitab25_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        MitabInputStream mitabInput = new MitabInputStream(solrResults, PsicquicSolrServer.DATA_FIELDS_25);

        Collection<psidev.psi.mi.tab.model.BinaryInteraction> binaryInteractions = mitabReader.read(mitabInput);
        Assert.assertNotNull(binaryInteractions);
        Assert.assertEquals(4, binaryInteractions.size());
    }

    @Test
    public void test_create_mitab26_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        MitabInputStream mitabInput = new MitabInputStream(solrResults, PsicquicSolrServer.DATA_FIELDS_26);

        Collection<psidev.psi.mi.tab.model.BinaryInteraction> binaryInteractions = mitabReader.read(mitabInput);
        Assert.assertNotNull(binaryInteractions);
        Assert.assertEquals(4, binaryInteractions.size());
    }

    @Test
    public void test_create_mitab27_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        MitabInputStream mitabInput = new MitabInputStream(solrResults, PsicquicSolrServer.DATA_FIELDS_27);

        Collection<psidev.psi.mi.tab.model.BinaryInteraction> binaryInteractions = mitabReader.read(mitabInput);
        Assert.assertNotNull(binaryInteractions);
        Assert.assertEquals(4, binaryInteractions.size());
    }
}
