package org.hupo.psi.mi.psicquic.indexing.batch.writer;

import org.apache.solr.client.solrj.SolrQuery;
import org.hupo.psi.calimocho.model.Row;
import org.hupo.psi.mi.psicquic.indexing.batch.AbstractSolrServerTest;
import org.hupo.psi.mi.psicquic.indexing.batch.reader.MitabCalimochoLineMapper;
import org.junit.Assert;
import org.junit.Test;
import org.springframework.batch.item.ExecutionContext;

import java.util.Arrays;

/**
 * Unit tester for SolrItemWriter
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>16/07/12</pre>
 */

public class SolrItemWriterUnitTest extends AbstractSolrServerTest {

    @Test
    public void test_write_rows() throws Exception {
        String solrUr= "http://localhost:18080/solr";
        SolrItemWriter writer = new SolrItemWriter();
        writer.setSolrUrl(solrUr);

        // add some data to the solrServer using writer
        MitabCalimochoLineMapper mitabLineMapper = new MitabCalimochoLineMapper();

        String mitab27 = "uniprotkb:P73045\tintact:EBI-1579103\tintact:EBI-1607518\tintact:EBI-1607516\tuniprotkb:slr1767(locus name)\tuniprotkb:alias2(gene name)\tpsi-mi:\"MI:0018\"(two hybrid)\tauthor et al.(2007)\tpubmed:18000013\ttaxid:4932(yeasx)|taxid:4932(\"Saccharomyces cerevisiae (Baker's yeast)\")\ttaxid:1142(9sync)|taxid:1142(Synechocystis)\tpsi-mi:\"MI:0915\"(physical association)\tpsi-mi:\"MI:0469\"(IntAct)\tintact:EBI-1607514\tauthor-score:C\t-\tpsi-mi:\"MI:0499\"(unspecified role)\tpsi-mi:\"MI:0499\"(unspecified role)\tpsi-mi:\"MI:0496\"(bait)\tpsi-mi:\"MI:0498\"(prey)\tpsi-mi:\"MI:0326\"(protein)\tpsi-mi:\"MI:0326\"(protein)\trefseq:NP_440386.1\t-\t-\t-\t-\t-\t-\t-\t2008/01/14\t2008/09/22\tcrc64:9E0E98F314F90177\t-\tintact-crc:3E26AC3853066993\t-\t-\t-\t-\t-\tpsi-mi:\"MI:0078\"(nucleotide sequence identification)\tpsi-mi:\"MI:0078\"(nucleotide sequence identification)\n";
        Row row = mitabLineMapper.mapLine(mitab27, 0);

        // index data to be hosted by PSICQUIC : we should have one result
        writer.open(new ExecutionContext());
        writer.write(Arrays.asList(row));
        writer.close();

        Assert.assertEquals(1L, solrJettyRunner.getSolrServer().query(new SolrQuery("*:*")).getResults().getNumFound());
    }
}
