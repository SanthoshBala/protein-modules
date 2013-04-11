package org.hupo.psi.mi.psicquic.model;

import org.apache.solr.client.solrj.SolrServer;
import org.apache.solr.client.solrj.SolrServerException;
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
import psidev.psi.mi.xml.converter.ConverterException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Unit tester of PSICQUICSolrServer
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/07/12</pre>
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration( locations = {"classpath*:/META-INF/psicquic-spring.xml",
        "classpath*:/jobs/psicquic-indexing-spring-test.xml"})
public class PsicquicSolrServerTest extends AbstractSolrServerTest {

    @Autowired
    private SolrMitabIndexer solrMitabIndexer;

    @Test
    public void test_wild_query() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        PsicquicSolrServer psicquicServer = new PsicquicSolrServer(server);

        PsicquicSearchResults results = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results);
        Assert.assertEquals(2L, results.getNumberResults());
    }

    @Test
    public void test_free_text_query_default_fields() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        PsicquicSolrServer psicquicServer = new PsicquicSolrServer(server);

        // default fields idA and idB
        PsicquicSearchResults results = psicquicServer.search("P07228", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results);
        Assert.assertEquals(1L, results.getNumberResults());

        // default fields altidA and altidB
        PsicquicSearchResults results2 = psicquicServer.search("EBI-5606437", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results2);
        Assert.assertEquals(1L, results2.getNumberResults());

        // default fields aliasA and aliasB
        PsicquicSearchResults results3 = psicquicServer.search("ITGB1", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results3);
        Assert.assertEquals(2L, results3.getNumberResults());

        // default fields pub id
        PsicquicSearchResults results4 = psicquicServer.search("9722563", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results4);
        Assert.assertEquals(2L, results4.getNumberResults());

        // default fields pub author
        PsicquicSearchResults results51 = psicquicServer.search("\"Loo DT et al.(1998)\" AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results51);
        Assert.assertEquals(1L, results51.getNumberResults());
        /*PsicquicSearchResults results52 = psicquicServer.search("\"Loo DT et al.\"", null, null, null, null);

        Assert.assertNotNull(results52);
        Assert.assertEquals(2L, results52.getNumberResults());
        PsicquicSearchResults results53 = psicquicServer.search("1998", null, null, null, null);

        Assert.assertNotNull(results53);
        Assert.assertEquals(1L, results53.getNumberResults());*/

        // default fields species
        PsicquicSearchResults results61 = psicquicServer.search("human-jurkat", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results61);
        Assert.assertEquals(1L, results61.getNumberResults());
        PsicquicSearchResults results62 = psicquicServer.search("9606", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results62);
        Assert.assertEquals(2L, results62.getNumberResults());

        // default fields detection method
        PsicquicSearchResults results7 = psicquicServer.search("\"two hybrid\" AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results7);
        Assert.assertEquals(1L, results7.getNumberResults());

        // default fields interaction type
        PsicquicSearchResults results8 = psicquicServer.search("\"direct interaction\" AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results8);
        Assert.assertEquals(1L, results8.getNumberResults());

        // default fields interaction identifier
        PsicquicSearchResults results9 = psicquicServer.search("EBI-5630468 AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results9);
        Assert.assertEquals(1L, results9.getNumberResults());

        // free text search, at least one which matches
        PsicquicSearchResults results10 = psicquicServer.search("P21333 P05556", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results10);
        Assert.assertEquals(2L, results10.getNumberResults());
    }

    @Test
    public void query_subset_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        PsicquicSolrServer psicquicServer = new PsicquicSolrServer(server);
        PsimiTabReader mitabReader = new PsimiTabReader();

        // query start is 2. Only 2 results are expected
        PsicquicSearchResults results = psicquicServer.search("negative:(true OR false)", 2, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results);
        Assert.assertEquals(2L, mitabReader.read(results.getMitab()).size());

        // query max is 3. Only 3 results are expected
        PsicquicSearchResults results2 = psicquicServer.search("negative:(true OR false)", null, 3, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results2);
        Assert.assertEquals(3L, mitabReader.read(results2.getMitab()).size());

        // query rows is 1 and query start = 2. Only 1 results are expected
        PsicquicSearchResults results3 = psicquicServer.search("negative:(true OR false)", 2, 1, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results3);
        Assert.assertEquals(1L, mitabReader.read(results3.getMitab()).size());

        // query rows is 2 and query start = 3. Only 1 results are expected because number max of results = 4
        PsicquicSearchResults results4 = psicquicServer.search("negative:(true OR false)", 3, 2, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results4);
        Assert.assertEquals(1L, mitabReader.read(results4.getMitab()).size());
    }

    @Test
    public void query_using_filters() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        PsicquicSolrServer psicquicServer = new PsicquicSolrServer(server);
        PsimiTabReader mitabReader = new PsimiTabReader();

        // identifier:P07228 is matching two results but with the query filter on participant detection method, only one result is expected
        PsicquicSearchResults results = psicquicServer.search("identifier:P07228 AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, "pmethod:\"predetermined participant\"");

        Assert.assertNotNull(results);
        Assert.assertEquals(1, mitabReader.read(results.getMitab()).size());

        // identifier:P21333 is matching four results but with the query filter on detection method + on negative, only one result is expected
        PsicquicSearchResults results2 = psicquicServer.searchWithFilters("identifier:P21333 AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, new String[]{"detmethod:\"anti bait coimmunoprecipitation\"", "negative:true"});

        Assert.assertNotNull(results2);
        Assert.assertEquals(1, mitabReader.read(results2.getMitab()).size());

        // 9722563 is matching four results (default query in pubid is working) but with the query filter on detection method + on negative, only one result is expected
        PsicquicSearchResults results3 = psicquicServer.searchWithFilters("9722563 AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, new String[]{"detmethod:\"anti bait coimmunoprecipitation\"", "negative:true"});

        Assert.assertNotNull(results3);
        Assert.assertEquals(1, mitabReader.read(results3.getMitab()).size());

        // western blot is not matching two results (default query in pmethod is not allowed)
        PsicquicSearchResults results4 = psicquicServer.searchWithFilters("\"western blot\" AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, new String[]{"detmethod:\"anti bait coimmunoprecipitation\"", "negative:true"});

        Assert.assertNotNull(results4);
        Assert.assertEquals(0, mitabReader.read(results4.getMitab()).size());
    }

    @Test
    public void test_valid_return_type() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, ConverterException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        PsicquicSolrServer psicquicServer = new PsicquicSolrServer(server);

        // default return type is mitab 25
        PsicquicSearchResults results = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB25, null);

        Assert.assertNotNull(results);
        InputStream mitab = results.getMitab();
        BufferedReader reader = new BufferedReader(new InputStreamReader(mitab));
        String firstLine = reader.readLine();
        reader.close();
        Assert.assertEquals(15, firstLine.split("\t").length);

        // return type is mitab 25
        PsicquicSearchResults results2 = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB25, null);

        Assert.assertNotNull(results2);
        InputStream mitab2 = results2.getMitab();
        BufferedReader reader2 = new BufferedReader(new InputStreamReader(mitab2));
        String firstLine2 = reader2.readLine();
        reader2.close();
        Assert.assertEquals(15, firstLine2.split("\t").length);

        // return type is mitab 26
        PsicquicSearchResults results3 = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB26, null);

        Assert.assertNotNull(results3);
        InputStream mitab3 = results3.getMitab();
        BufferedReader reader3 = new BufferedReader(new InputStreamReader(mitab3));
        String firstLine3 = reader3.readLine();
        reader3.close();
        Assert.assertEquals(36, firstLine3.split("\t").length);

        // return type is mitab 27
        PsicquicSearchResults results4 = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results4);
        InputStream mitab4 = results4.getMitab();
        BufferedReader reader4 = new BufferedReader(new InputStreamReader(mitab4));
        String firstLine4 = reader4.readLine();
        reader4.close();
        Assert.assertEquals(42, firstLine4.split("\t").length);

        // return type is mitab 27 for return type xml
        PsicquicSearchResults results5 = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_XML25, null);

        Assert.assertNotNull(results5);
        InputStream mitab5 = results5.getMitab();
        BufferedReader reader5 = new BufferedReader(new InputStreamReader(mitab5));
        String firstLine5 = reader5.readLine();
        reader5.close();
        Assert.assertEquals(42, firstLine5.split("\t").length);

        // return type is count, no mitab is generated
        PsicquicSearchResults results6 = psicquicServer.search("*", null, null, PsicquicSolrServer.RETURN_TYPE_COUNT, null);

        Assert.assertNotNull(results6);
        Assert.assertNull(results6.getMitab());
    }

    @Test
    public void test_wild_query_to_lower_case() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        PsicquicSolrServer psicquicServer = new PsicquicSolrServer(server);

        // should return 4 lines because P0* can match both P07228 and P05556
        PsicquicSearchResults results = psicquicServer.search("P0* AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results);
        Assert.assertEquals(4L, results.getNumberResults());

        // should return 4 lines because P0* can match both P07228 and P05556
        PsicquicSearchResults results2 = psicquicServer.search("identifier:P0* AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results2);
        Assert.assertEquals(4L, results2.getNumberResults());

        // should return 2 lines because P0* can match both P07228 and P05556
        PsicquicSearchResults results3 = psicquicServer.search("idA:P0* AND negative:(true OR false)", null, null, PsicquicSolrServer.RETURN_TYPE_MITAB27, null);

        Assert.assertNotNull(results3);
        Assert.assertEquals(2L, results3.getNumberResults());
    }
}
