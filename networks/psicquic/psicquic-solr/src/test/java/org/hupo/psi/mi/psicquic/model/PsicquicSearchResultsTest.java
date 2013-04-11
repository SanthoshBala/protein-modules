package org.hupo.psi.mi.psicquic.model;

import org.apache.solr.client.solrj.SolrQuery;
import org.apache.solr.client.solrj.SolrServer;
import org.apache.solr.client.solrj.SolrServerException;
import org.apache.solr.common.SolrDocumentList;
import org.hupo.psi.calimocho.io.IllegalRowException;
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
import psidev.psi.mi.tab.converter.tab2xml.XmlConversionException;
import psidev.psi.mi.tab.model.BinaryInteraction;
import psidev.psi.mi.xml.model.Entry;
import psidev.psi.mi.xml.model.EntrySet;

import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;

/**
 * Unit tester of PsicquicSearchResult
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/07/12</pre>
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration( locations = {"classpath*:/META-INF/psicquic-spring.xml",
        "classpath*:/jobs/psicquic-indexing-spring-test.xml"})
public class PsicquicSearchResultsTest  extends AbstractSolrServerTest {

    @Autowired
    private SolrMitabIndexer solrMitabIndexer;

    @Test
    public void test_create_mitab_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        PsicquicSearchResults psicquicSearchResults = new PsicquicSearchResults(solrResults, PsicquicSolrServer.DATA_FIELDS_27);

        InputStream mitabInputStream = psicquicSearchResults.getMitab();
        Assert.assertNotNull(mitabInputStream);

        PsimiTabReader mitabReader = new PsimiTabReader();
        Collection<BinaryInteraction> binaryInteractions = mitabReader.read(mitabInputStream);
        Assert.assertNotNull(binaryInteractions);
        Assert.assertEquals(4, binaryInteractions.size());

        mitabInputStream.close();

        PsicquicSearchResults psicquicSearchResults2 = new PsicquicSearchResults(solrResults, new String[] {});

        InputStream mitabInputStream2 = psicquicSearchResults2.getMitab();
        Assert.assertNull(mitabInputStream2);
    }

    @Test
    public void test_create_xml_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException, XmlConversionException, IllegalAccessException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        PsicquicSearchResults psicquicSearchResults = new PsicquicSearchResults(solrResults, PsicquicSolrServer.DATA_FIELDS_27);

        EntrySet xml = psicquicSearchResults.createEntrySet();
        Assert.assertNotNull(xml);

        Collection<Entry> entries = xml.getEntries();
        Assert.assertEquals(1, entries.size());
        Entry entry = entries.iterator().next();
        Assert.assertEquals(4, entry.getInteractions().size());
    }

    @Test
    public void test_create_xgmml_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException, XmlConversionException, IllegalAccessException, IllegalRowException, XMLStreamException, JAXBException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        PsicquicSearchResults psicquicSearchResults = new PsicquicSearchResults(solrResults, PsicquicSolrServer.DATA_FIELDS_27);

        InputStream xgmml = psicquicSearchResults.createXGMML();
        Assert.assertNotNull(xgmml);

        xgmml.close();
    }

    @Test
    public void test_create_rdf_results() throws JobInstanceAlreadyCompleteException, JobParametersInvalidException, JobRestartException, JobExecutionAlreadyRunningException, SolrServerException, PsicquicSolrException, PsimiTabException, IOException, XmlConversionException, IllegalAccessException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList solrResults = server.query(new SolrQuery("*:*")).getResults();
        PsicquicSearchResults psicquicSearchResults = new PsicquicSearchResults(solrResults, PsicquicSolrServer.DATA_FIELDS_27);

        InputStream invalidFormat = psicquicSearchResults.createRDFOrBiopax("rdf");
        Assert.assertNull(invalidFormat);

        InputStream biopaxFormat = psicquicSearchResults.createRDFOrBiopax("biopax");
        Assert.assertNotNull(biopaxFormat);

        biopaxFormat.close();

        InputStream biopaxL3Format = psicquicSearchResults.createRDFOrBiopax("biopax-L3");
        Assert.assertNotNull(biopaxL3Format);

        biopaxL3Format.close();

        InputStream biopaxL2Format = psicquicSearchResults.createRDFOrBiopax("biopax-L2");
        Assert.assertNotNull(biopaxL2Format);

        biopaxL2Format.close();

        InputStream rdfXmlFormat = psicquicSearchResults.createRDFOrBiopax("rdf-xml");
        Assert.assertNotNull(rdfXmlFormat);

        rdfXmlFormat.close();

        InputStream rdfXmlAbbrevFormat = psicquicSearchResults.createRDFOrBiopax("rdf-xml-abbrev");
        Assert.assertNotNull(rdfXmlAbbrevFormat);

        rdfXmlAbbrevFormat.close();

        InputStream rdfN3 = psicquicSearchResults.createRDFOrBiopax("rdf-n3");
        Assert.assertNotNull(rdfN3);

        rdfN3.close();

        InputStream rdfTurtle = psicquicSearchResults.createRDFOrBiopax("rdf-turtle");
        Assert.assertNotNull(rdfTurtle);

        rdfTurtle.close();
    }
}
