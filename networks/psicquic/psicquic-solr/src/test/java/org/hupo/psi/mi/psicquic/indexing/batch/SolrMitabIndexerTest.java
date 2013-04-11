package org.hupo.psi.mi.psicquic.indexing.batch;

import org.apache.solr.client.solrj.SolrQuery;
import org.apache.solr.client.solrj.SolrServer;
import org.apache.solr.client.solrj.impl.HttpSolrServer;
import org.apache.solr.common.SolrDocument;
import org.apache.solr.common.SolrDocumentList;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.batch.core.repository.JobExecutionAlreadyRunningException;
import org.springframework.batch.core.repository.JobRestartException;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationContext;
import org.springframework.test.annotation.DirtiesContext;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;
import psidev.psi.mi.calimocho.solr.converter.SolrFieldName;

import java.util.Collection;
import java.util.Iterator;

/**
 * TODO comment this
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/05/12</pre>
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration( locations = {"classpath*:/META-INF/psicquic-spring.xml",
        "classpath*:/jobs/psicquic-indexing-spring-test.xml"})
public class SolrMitabIndexerTest extends AbstractSolrServerTest {

    @Autowired
    private SolrMitabIndexer solrMitabIndexer;

    @Autowired
    private ApplicationContext applicationContext;

    @Test
    public void test_indexing_negative() throws Exception, JobRestartException, JobExecutionAlreadyRunningException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        HttpSolrServer server = solrJettyRunner.getSolrServer();

        Assert.assertEquals(4L, server.query(new SolrQuery("*:*")).getResults().getNumFound());
        // two negative, two positive
        Assert.assertEquals(2L, server.query(new SolrQuery("negative:true")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("negative:false")).getResults().getNumFound());

        // test idA and idB -> identifier and id
        Assert.assertEquals(2L, server.query(new SolrQuery("id:P07228")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("identifier:P07228")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("idA:P07228")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("idB:P07228")).getResults().getNumFound());

        // test altidA and altidB -> identifier and id
        Assert.assertEquals(2L, server.query(new SolrQuery("id:EBI-5606437")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("identifier:EBI-5606437")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idA:EBI-5606437")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idB:EBI-5606437")).getResults().getNumFound());

        // test altiasA and aliasB -> identifier and alias
        Assert.assertEquals(2L, server.query(new SolrQuery("alias:RGD-receptor")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("identifier:RGD-receptor")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idA:RGD-receptor")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idB:RGD-receptor")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("id:RGD-receptor")).getResults().getNumFound());

        // test detmethod
        Assert.assertEquals(1L, server.query(new SolrQuery("detmethod:\"two hybrid\"")).getResults().getNumFound());

        // test author
        Assert.assertEquals(1L, server.query(new SolrQuery("pubauth:\"Loo DT et al.(1998)\"")).getResults().getNumFound());

        // test pub id
        Assert.assertEquals(4L, server.query(new SolrQuery("pubid:9722563")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pubid:imex")).getResults().getNumFound());

        // test protein species
        Assert.assertEquals(2L, server.query(new SolrQuery("taxidA:\"human-jurkat\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("taxidB:\"human-jurkat\"")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("species:\"human-jurkat\"")).getResults().getNumFound());

        // test interaction type
        Assert.assertEquals(1L, server.query(new SolrQuery("type:\"MI:0407\"")).getResults().getNumFound());

        // test interaction id
        Assert.assertEquals(1L, server.query(new SolrQuery("interaction_id:EBI-5630468")).getResults().getNumFound());

        // test biorole
        Assert.assertEquals(4L, server.query(new SolrQuery("pbioroleA:\"unspecified role\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbioroleB:\"unspecified role\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbiorole:\"unspecified role\"")).getResults().getNumFound());

        // test interactor type
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeA:protein")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeB:protein")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptype:protein")).getResults().getNumFound());

        // test interactor type
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeA:protein")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeB:protein")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptype:protein")).getResults().getNumFound());

        // test interactor xref
        Assert.assertEquals(1L, server.query(new SolrQuery("pxrefA:\"GO:0008305\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pxrefB:\"GO:0008305\"")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("pxref:\"GO:0008305\"")).getResults().getNumFound());

        // test interaction xref
        Assert.assertEquals(1L, server.query(new SolrQuery("xref:\"IM-17229-4\"")).getResults().getNumFound());

        // test interaction annotations
        Assert.assertEquals(1L, server.query(new SolrQuery("annot:\"Fig. 1.\"")).getResults().getNumFound());

        // test udate
        Assert.assertEquals(1L, server.query(new SolrQuery("udate:[20120301 TO 20120302]")).getResults().getNumFound());

        // test feature
        Assert.assertEquals(1L, server.query(new SolrQuery("ftypeA:\"necessary binding region\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("ftypeB:\"necessary binding region\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("ftype:\"necessary binding region\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("ftypeA:2171-2647")).getResults().getNumFound());

        // test stc
        Assert.assertEquals(0L, server.query(new SolrQuery("stc:true")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("stc:false")).getResults().getNumFound());

        // test participant detection method
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethodA:\"MI:0981\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethodB:\"MI:0981\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethod:\"MI:0981\"")).getResults().getNumFound());
    }

    @Test
    public void test_indexing_parameter() throws Exception, JobRestartException, JobExecutionAlreadyRunningException {

        solrMitabIndexer.startJob("mitabIndexParameterJob");

        HttpSolrServer server = solrJettyRunner.getSolrServer();

        Assert.assertEquals(14L, server.query(new SolrQuery("*:*")).getResults().getNumFound());
        // two parameters
        Assert.assertEquals(3L, server.query(new SolrQuery("param:true")).getResults().getNumFound());
        Assert.assertEquals(11L, server.query(new SolrQuery("param:false")).getResults().getNumFound());

        // test spoke expansion
        Assert.assertEquals(8L, server.query(new SolrQuery("complex:\"MI:1060\"")).getResults().getNumFound());
        Assert.assertEquals(8L, server.query(new SolrQuery("complex:\"spoke expansion\"")).getResults().getNumFound());
        Assert.assertEquals(8L, server.query(new SolrQuery("complex:\"psi-mi:MI:1060\"")).getResults().getNumFound());
        //psi-mi:"MI:1060"(spoke expansion)

        // test stc
        Assert.assertEquals(3L, server.query(new SolrQuery("stc:true")).getResults().getNumFound());
        Assert.assertEquals(11L, server.query(new SolrQuery("stc:false")).getResults().getNumFound());
    }

    @Test
    public void test_search_different_values() throws Exception, JobRestartException, JobExecutionAlreadyRunningException {

        solrMitabIndexer.startJob("mitabIndexNegativeJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        Assert.assertEquals(4L, server.query(new SolrQuery("*:*")).getResults().getNumFound());

        // test idA and idB -> identifier and id
        Assert.assertEquals(4L, server.query(new SolrQuery("id:uniprotkb")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("identifier:uniprotkb")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("idA:uniprotkb")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("idB:uniprotkb")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("id:\"uniprotkb:P07228\"")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("identifier:\"uniprotkb:P07228\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("idA:\"uniprotkb:P07228\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("idB:\"uniprotkb:P07228\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idB:\"ensembl:P07228\"")).getResults().getNumFound());

        // test altidA and altidB -> identifier and id
        Assert.assertEquals(4L, server.query(new SolrQuery("id:intact")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("identifier:intact")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("id:\"intact:EBI-5606437\"")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("identifier:\"intact:EBI-5606437\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("identifier:\"uniprotkb:EBI-5606437\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idA:\"intact:EBI-5606437\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idB:\"intact:EBI-5606437\"")).getResults().getNumFound());

        // test altiasA and aliasB -> identifier and alias
        Assert.assertEquals(4L, server.query(new SolrQuery("alias:uniprotkb")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("alias:\"uniprotkb:RGD-receptor\"")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("identifier:\"uniprotkb:RGD-receptor\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("alias:\"intact:RGD-receptor\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idA:\"uniprotkb:RGD-receptor\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("idB:\"uniprotkb:RGD-receptor\"")).getResults().getNumFound());
        Assert.assertEquals(0L, server.query(new SolrQuery("id:\"uniprotkb:RGD-receptor\"")).getResults().getNumFound());

        // test detmethod
        Assert.assertEquals(4L, server.query(new SolrQuery("detmethod:\"psi-mi\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("detmethod:\"MI:0018\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("detmethod:\"psi-mi:MI:0018\"")).getResults().getNumFound());

        // test pub id
        Assert.assertEquals(4L, server.query(new SolrQuery("pubid:\"pubmed:9722563\"")).getResults().getNumFound());

        // test protein species
        Assert.assertEquals(3L, server.query(new SolrQuery("taxidA:9606")).getResults().getNumFound());
        Assert.assertEquals(3L, server.query(new SolrQuery("taxidB:\"taxid:9606\"")).getResults().getNumFound());

        // test interaction type
        Assert.assertEquals(1L, server.query(new SolrQuery("type:\"direct interaction\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("type:\"psi-mi:MI:0407\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("type:psi-mi")).getResults().getNumFound());

        // test interaction id
        Assert.assertEquals(1L, server.query(new SolrQuery("interaction_id:\"intact:EBI-5630468\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("interaction_id:intact")).getResults().getNumFound());

        // test biorole
        Assert.assertEquals(4L, server.query(new SolrQuery("pbioroleA:\"psi-mi:MI:0499\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbioroleB:\"psi-mi:MI:0499\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbiorole:\"psi-mi:MI:0499\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbioroleA:\"psi-mi:MI:0499\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbioroleB:psi-mi")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pbiorole:psi-mi")).getResults().getNumFound());

        // test interactor type
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeA:psi-mi")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeB:psi-mi")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptype:psi-mi")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeA:\"psi-mi:MI:0326\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptypeB:\"psi-mi:MI:0326\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("ptype:\"psi-mi:MI:0326\"")).getResults().getNumFound());

        // test interactor xref
        Assert.assertEquals(4L, server.query(new SolrQuery("pxrefA:go")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pxrefB:go")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pxref:go")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pxrefA:\"go:GO:0008305\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pxrefB:\"go:GO:0008305\"")).getResults().getNumFound());
        Assert.assertEquals(2L, server.query(new SolrQuery("pxref:\"go:GO:0008305\"")).getResults().getNumFound());

        // test interaction xref
        Assert.assertEquals(1L, server.query(new SolrQuery("xref:\"imex:IM-17229-4\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("xref:imex")).getResults().getNumFound());

        // test interaction annotations
        Assert.assertEquals(4L, server.query(new SolrQuery("annot:\"figure legend\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("annot:\"figure legend:Fig. 1.\"")).getResults().getNumFound());

        // test udate
        Assert.assertEquals(1L, server.query(new SolrQuery("udate:20120301")).getResults().getNumFound());

        // test participant detection method
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethodA:\"psi-mi:MI:0981\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethodB:\"psi-mi:MI:0981\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethod:\"psi-mi:MI:0981\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethodA:\"tag visualisation by peroxidase activity\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethodB:\"tag visualisation by peroxidase activity\"")).getResults().getNumFound());
        Assert.assertEquals(1L, server.query(new SolrQuery("pmethod:\"tag visualisation by peroxidase activity\"")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pmethodA:psi-mi")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pmethodB:psi-mi")).getResults().getNumFound());
        Assert.assertEquals(4L, server.query(new SolrQuery("pmethod:psi-mi")).getResults().getNumFound());
    }

    @Test
    public void test_retrieve_results_solr() throws Exception, JobExecutionAlreadyRunningException {

        solrMitabIndexer.startJob("mitabIndexParameterJob");

        SolrServer server = solrJettyRunner.getSolrServer();

        SolrDocumentList results = server.query(new SolrQuery("interaction_id:EBI-4321313")).getResults();

        Iterator<SolrDocument> iterator = results.iterator();

        Assert.assertEquals(2, results.getNumFound());
        SolrDocument doc = iterator.next();

        Collection<Object> idA = doc.getFieldValues(SolrFieldName.idA.toString()+"_o");
        Assert.assertEquals(1, idA.size());
        Assert.assertEquals("uniprotkb:O43521", idA.iterator().next());
        Collection<Object> idB = doc.getFieldValues(SolrFieldName.idB.toString()+"_o");
        Assert.assertEquals(1, idB.size());
        Assert.assertEquals("uniprotkb:P10415", idB.iterator().next());
        Collection<Object> altidA = doc.getFieldValues(SolrFieldName.altidA.toString()+"_o");
        Assert.assertEquals(1, altidA.size());
        Assert.assertEquals("intact:EBI-526406", altidA.iterator().next());
        Collection<Object> altidB = doc.getFieldValues(SolrFieldName.altidB.toString()+"_o");
        Assert.assertEquals(1, altidB.size());
        Assert.assertEquals("intact:EBI-77694", altidB.iterator().next());
        Collection<Object> aliasA = doc.getFieldValues(SolrFieldName.aliasA.toString()+"_o");
        Assert.assertEquals(1, aliasA.size());
        Assert.assertEquals("uniprotkb:BCL2L11(gene name)|uniprotkb:BIM(gene name synonym)|uniprotkb:Bcl2-interacting mediator of cell death(gene name synonym)", aliasA.iterator().next());
        Collection<Object> aliasB = doc.getFieldValues(SolrFieldName.aliasB.toString()+"_o");
        Assert.assertEquals(1, aliasB.size());
        Assert.assertEquals("uniprotkb:BCL2(gene name)", aliasB.iterator().next());
        Collection<Object> detmethod = doc.getFieldValues(SolrFieldName.detmethod.toString()+"_o");
        Assert.assertEquals(1, detmethod.size());
        Assert.assertEquals("psi-mi:\"MI:0107\"(surface plasmon resonance)", detmethod.iterator().next());
        Collection<Object> pubauth = doc.getFieldValues(SolrFieldName.pubauth.toString()+"_o");
        Assert.assertEquals(1, pubauth.size());
        Assert.assertEquals("author et al.(2008)", pubauth.iterator().next());
        Collection<Object> pubid = doc.getFieldValues(SolrFieldName.pubid.toString()+"_o");
        Assert.assertEquals(1, pubid.size());
        Assert.assertEquals("pubmed:21199865|imex:IM-16869", pubid.iterator().next());
        Collection<Object> taxidA = doc.getFieldValues(SolrFieldName.taxidA.toString()+"_o");
        Assert.assertEquals(1, taxidA.size());
        Assert.assertEquals("taxid:-1(in vitro)|taxid:-1(\"In vitro (In vitro)\")|taxid:9606(human)|taxid:9606(Homo sapiens)", taxidA.iterator().next());
        Collection<Object> taxidB = doc.getFieldValues(SolrFieldName.taxidB.toString()+"_o");
        Assert.assertEquals(1, taxidB.size());
        Assert.assertEquals("taxid:9606(human)|taxid:9606(Homo sapiens)", taxidB.iterator().next());
        Collection<Object> types = doc.getFieldValues(SolrFieldName.type.toString()+"_o");
        Assert.assertEquals(1, types.size());
        Assert.assertEquals("psi-mi:\"MI:0915\"(physical association)", types.iterator().next());
        Collection<Object> source = doc.getFieldValues(SolrFieldName.source.toString()+"_o");
        Assert.assertEquals(1, source.size());
        Assert.assertEquals("psi-mi:\"MI:0469\"(IntAct)", source.iterator().next());
        Collection<Object> interaction_id = doc.getFieldValues(SolrFieldName.interaction_id.toString()+"_o");
        Assert.assertEquals(1, interaction_id.size());
        Assert.assertEquals("intact:EBI-4321313", interaction_id.iterator().next());
        Collection<Object> confidence = doc.getFieldValues(SolrFieldName.confidence.toString()+"_o");
        Assert.assertEquals(1, confidence.size());
        Assert.assertEquals("author-score:high", confidence.iterator().next());
        Collection<Object> complex = doc.getFieldValues(SolrFieldName.complex.toString()+"_o");
        Assert.assertEquals(1, complex.size());
        Assert.assertEquals("psi-mi:\"MI:1060\"(spoke expansion)", complex.iterator().next());
        Collection<Object> bioRoleA = doc.getFieldValues(SolrFieldName.pbioroleA.toString()+"_o");
        Assert.assertEquals(1, bioRoleA.size());
        Assert.assertEquals("psi-mi:\"MI:0499\"(unspecified role)", bioRoleA.iterator().next());
        Collection<Object> bioRoleB = doc.getFieldValues(SolrFieldName.pbioroleB.toString()+"_o");
        Assert.assertEquals(1, bioRoleB.size());
        Assert.assertEquals("psi-mi:\"MI:0499\"(unspecified role)", bioRoleB.iterator().next());
        Collection<Object> expRoleA = doc.getFieldValues(SolrFieldName.pexproleA.toString()+"_o");
        Assert.assertEquals(1, expRoleA.size());
        Assert.assertEquals("psi-mi:\"MI:0496\"(bait)", expRoleA.iterator().next());
        Collection<Object> expRoleB = doc.getFieldValues(SolrFieldName.pexproleB.toString()+"_o");
        Assert.assertEquals(1, expRoleB.size());
        Assert.assertEquals("psi-mi:\"MI:0498\"(prey)", expRoleB.iterator().next());
        Collection<Object> ptypeA = doc.getFieldValues(SolrFieldName.ptypeA.toString()+"_o");
        Assert.assertEquals(1, ptypeA.size());
        Assert.assertEquals("psi-mi:\"MI:0326\"(protein)", ptypeA.iterator().next());
        Collection<Object> ptypeB = doc.getFieldValues(SolrFieldName.ptypeB.toString()+"_o");
        Assert.assertEquals(1, ptypeB.size());
        Assert.assertEquals("psi-mi:\"MI:0326\"(protein)", ptypeB.iterator().next());
        Collection<Object> pxrefA = doc.getFieldValues(SolrFieldName.pxrefA.toString()+"_o");
        Assert.assertEquals(1, pxrefA.size());
        Assert.assertEquals("uniprotkb:Q0MSE7(secondary-ac)|uniprotkb:Q0MSE8(secondary-ac)|uniprotkb:Q0MSE9(secondary-ac)|uniprotkb:Q53R28(secondary-ac)|uniprotkb:Q6JTU6(secondary-ac)|uniprotkb:Q6T851(secondary-ac)|uniprotkb:Q6TE14(secondary-ac)|uniprotkb:Q6TE15(secondary-ac)|uniprotkb:Q6TE16(secondary-ac)|uniprotkb:Q6V402(secondary-ac)|uniprotkb:A8K2W2(secondary-ac)|uniprotkb:Q8WYL6(secondary-ac)|uniprotkb:Q8WYL7(secondary-ac)|uniprotkb:Q8WYL8(secondary-ac)|uniprotkb:Q8WYL9(secondary-ac)|uniprotkb:Q8WYM0(secondary-ac)|uniprotkb:Q8WYM1(secondary-ac)|uniprotkb:O43522(secondary-ac)|intact:EBI-1002160(intact-secondary)|refseq:NP_001191035.1|refseq:NP_001191036.1|refseq:NP_001191037.1|refseq:NP_001191038.1|refseq:NP_001191039.1|refseq:NP_001191040.1|refseq:NP_001191041.1|refseq:NP_001191042.1|refseq:NP_006529.1|refseq:NP_619527.1|refseq:NP_619528.1|refseq:NP_619529.1|refseq:NP_619530.1|refseq:NP_619531.1|refseq:NP_619532.1|refseq:NP_619533.1|refseq:NP_996885.1|refseq:NP_996886.1|interpro:IPR014771|interpro:IPR017288|interpro:IPR015040|go:\"GO:0005829\"|go:\"GO:0012505\"|go:\"GO:0005741\"|go:\"GO:0005886\"|go:\"GO:0005515\"|go:\"GO:0008633\"|go:\"GO:0008624\"|go:\"GO:0008629\"|go:\"GO:0048011\"|go:\"GO:0032464\"|go:\"GO:0090200\"|ensembl:ENSG00000153094|ensembl:ENSG00000153094|reactome:REACT_578|ipi:IPI00012853|ipi:IPI00103743|ipi:IPI00216585|ipi:IPI00410159|ipi:IPI00428840|ipi:IPI00451139|ipi:IPI00514401|ipi:IPI00873921|ipi:IPI00878323|ipi:IPI00914559|ipi:IPI00914560|ipi:IPI00914563|ipi:IPI00914586|ipi:IPI00914592|ipi:IPI00914600|ipi:IPI00914670|ipi:IPI00914677|rcsb pdb:2K7W|rcsb pdb:2NL9|rcsb pdb:2V6Q|rcsb pdb:2VM6|rcsb pdb:2WH6|rcsb pdb:3D7V|rcsb pdb:3FDL|rcsb pdb:3IO8|rcsb pdb:3IO9|rcsb pdb:3KJ0|rcsb pdb:3KJ1|rcsb pdb:3KJ2|go:\"GO:2001244\"|go:\"GO:0097141\"|go:\"GO:0097140\"|go:\"GO:0006919\"|go:\"GO:0043065\"|reactome:REACT_111102", pxrefA.iterator().next());
        Collection<Object> pxrefB = doc.getFieldValues(SolrFieldName.pxrefB.toString()+"_o");
        Assert.assertEquals(1, pxrefB.size());
        Assert.assertEquals("go:\"GO:0001836\"|go:\"GO:0034097\"|go:\"GO:0006974\"|go:\"GO:0042493\"|go:\"GO:0010039\"|go:\"GO:0035094\"|go:\"GO:0009636\"|ensembl:ENSG00000171791|ensembl:ENSG00000171791|reactome:REACT_578|reactome:REACT_6900|ipi:IPI00020961|ipi:IPI00217817|rcsb pdb:1G5M|rcsb pdb:1GJH|rcsb pdb:1YSW|rcsb pdb:2O21|rcsb pdb:2O22|rcsb pdb:2O2F|rcsb pdb:2W3L|rcsb pdb:2XA0|go:\"GO:0043497\"|go:\"GO:0043496\"|go:\"GO:0022898\"|uniprotkb:P10416(secondary-ac)|uniprotkb:Q13842(secondary-ac)|uniprotkb:Q16197(secondary-ac)|refseq:NP_000624.2|interpro:IPR013278|interpro:IPR002475|interpro:IPR000712|interpro:IPR020717|interpro:IPR020726|interpro:IPR020728|interpro:IPR003093|interpro:IPR020731|interpro:IPR004725|go:\"GO:0005789\"|go:\"GO:0005741\"|go:\"GO:0031965\"|go:\"GO:0046930\"|go:\"GO:0051434\"|go:\"GO:0015267\"|go:\"GO:0002020\"|go:\"GO:0046982\"|go:\"GO:0042803\"|go:\"GO:0043565\"|go:\"GO:0031625\"|go:\"GO:0008633\"|go:\"GO:0006916\"|go:\"GO:0070059\"|go:\"GO:0042100\"|go:\"GO:0050853\"|go:\"GO:0051607\"|go:\"GO:0007565\"|go:\"GO:0006959\"|go:\"GO:0008629\"|go:\"GO:0045087\"|go:\"GO:0032848\"|go:\"GO:0051902\"|go:\"GO:0043524\"|go:\"GO:0051402\"|go:\"GO:0030890\"|go:\"GO:0030307\"|go:\"GO:0000209\"|go:\"GO:0046902\"|go:\"GO:0051881\"|go:\"GO:0005515\"|go:\"GO:0035872\"|go:\"GO:2001234\"|go:\"GO:0051924\"", pxrefB.iterator().next());
        Collection<Object> xref = doc.getFieldValues(SolrFieldName.xref.toString()+"_o");
        Assert.assertEquals(1, xref.size());
        Assert.assertEquals("imex:IM-16869-5(imex-primary)", xref.iterator().next());
        Collection<Object> annotA = doc.getFieldValues(SolrFieldName.annotA.toString()+"_o");
        Assert.assertEquals(1, annotA.size());
        Assert.assertEquals("caution:test", annotA.iterator().next());
        Collection<Object> annotB = doc.getFieldValues(SolrFieldName.annotB.toString()+"_o");
        Assert.assertEquals(1, annotB.size());
        Assert.assertEquals("comment:test", annotB.iterator().next());
        Collection<Object> annot = doc.getFieldValues(SolrFieldName.annot.toString()+"_o");
        Assert.assertEquals(1, annot.size());
        Assert.assertEquals("curation depth:imex curation|full coverage:Only protein-protein interactions|figure legend:Table 2", annot.iterator().next());
        Collection<Object> host = doc.getFieldValues(SolrFieldName.taxidHost.toString()+"_o");
        Assert.assertEquals(1, host.size());
        Assert.assertEquals("taxid:-1(in vitro)", host.iterator().next());
        Collection<Object> parameter = doc.getFieldValues(SolrFieldName.param.toString()+"_o");
        Assert.assertEquals(1, parameter.size());
        Assert.assertEquals("kd:0.8", parameter.iterator().next());
        Collection<Object> created = doc.getFieldValues(SolrFieldName.cdate.toString()+"_o");
        Assert.assertEquals(1, created.size());
        Assert.assertEquals("2011/08/03", created.iterator().next());
        Collection<Object> updated = doc.getFieldValues(SolrFieldName.udate.toString()+"_o");
        Assert.assertEquals(1, updated.size());
        Assert.assertEquals("2011/08/04", updated.iterator().next());
        Collection<Object> checksumA = doc.getFieldValues(SolrFieldName.checksumA.toString()+"_o");
        Assert.assertEquals(1, checksumA.size());
        Assert.assertEquals("crc64:D75735E469CA6997", checksumA.iterator().next());
        Collection<Object> checksumB = doc.getFieldValues(SolrFieldName.checksumB.toString()+"_o");
        Assert.assertEquals(1, checksumB.size());
        Assert.assertEquals("crc64:3C49F2B714DC9CCB", checksumB.iterator().next());
        Collection<Object> checksum = doc.getFieldValues(SolrFieldName.checksumI.toString()+"_o");
        Assert.assertEquals(1, checksum.size());
        Assert.assertEquals("intact-crc:xxxxx", checksum.iterator().next());
        Collection<Object> negative = doc.getFieldValues(SolrFieldName.negative.toString()+"_o");
        Assert.assertNull(negative);
        Collection<Object> ftypeA = doc.getFieldValues(SolrFieldName.ftypeA.toString()+"_o");
        Assert.assertEquals(1, ftypeA.size());
        Assert.assertEquals("sufficient binding region:141-166", ftypeA.iterator().next());
        Collection<Object> ftypeB = doc.getFieldValues(SolrFieldName.ftypeB.toString()+"_o");
        Assert.assertNull(ftypeB);
        Collection<Object> stcA = doc.getFieldValues(SolrFieldName.stcA.toString()+"_o");
        Assert.assertEquals(1, stcA.size());
        Assert.assertEquals("1", stcA.iterator().next());
        Collection<Object> stcB = doc.getFieldValues(SolrFieldName.stcB.toString()+"_o");
        Assert.assertEquals(1, stcB.size());
        Assert.assertEquals("2", stcB.iterator().next());
        Collection<Object> pmethodA = doc.getFieldValues(SolrFieldName.pmethodA.toString()+"_o");
        Assert.assertEquals(1, pmethodA.size());
        Assert.assertEquals("psi-mi:\"MI:0396\"(predetermined participant)", pmethodA.iterator().next());
        Collection<Object> pmethodB = doc.getFieldValues(SolrFieldName.pmethodB.toString()+"_o");
        Assert.assertEquals(1, pmethodB.size());
        Assert.assertEquals("psi-mi:\"MI:0396\"(predetermined participant)", pmethodB.iterator().next());
    }

    @Test
    @DirtiesContext
    public void failingRelease1() throws Exception {

        // first time should fail after 2 readings
        solrMitabIndexer.startJob("mitabIndexFailingJob");
        SolrServer server = solrJettyRunner.getSolrServer();

        Assert.assertEquals(2L, server.query(new SolrQuery("*:*")).getResults().getNumFound());

        // it should resume and continue
        solrMitabIndexer.resumeJob("mitabIndexFailingJob");

        Assert.assertEquals(4L, server.query(new SolrQuery("*:*")).getResults().getNumFound());
    }
}


