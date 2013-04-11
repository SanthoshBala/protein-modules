package org.hupo.psi.mi.psicquic.indexing.batch.reader;

import junit.framework.Assert;
import org.hupo.psi.calimocho.model.Row;
import org.junit.Test;

/**
 * This class is testing MitabCalimochoLineMapper
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>12/07/12</pre>
 */
public class MitabCalimochoLineMapperTest {

    @Test
    public void read_mitab_27() throws Exception {

        MitabCalimochoLineMapper mitabLineMapper = new MitabCalimochoLineMapper();

        String mitab27 = "uniprotkb:P73045\tintact:EBI-1579103\tintact:EBI-1607518\tintact:EBI-1607516\tuniprotkb:slr1767(locus name)\tuniprotkb:alias2(gene name)\tpsi-mi:\"MI:0018\"(two hybrid)\tauthor et al.(2007)\tpubmed:18000013\ttaxid:4932(yeasx)|taxid:4932(\"Saccharomyces cerevisiae (Baker's yeast)\")\ttaxid:1142(9sync)|taxid:1142(Synechocystis)\tpsi-mi:\"MI:0915\"(physical association)\tpsi-mi:\"MI:0469\"(IntAct)\tintact:EBI-1607514\tauthor-score:C\t-\tpsi-mi:\"MI:0499\"(unspecified role)\tpsi-mi:\"MI:0499\"(unspecified role)\tpsi-mi:\"MI:0496\"(bait)\tpsi-mi:\"MI:0498\"(prey)\tpsi-mi:\"MI:0326\"(protein)\tpsi-mi:\"MI:0326\"(protein)\trefseq:NP_440386.1\trefseq:test2\tirefindex:test3\tcomment:test1\tcomment:test2\tcomment:test3\ttaxid:9606(human)\tkd:1.5\t2008/01/14\t2008/09/22\tcrc64:9E0E98F314F90177\tcrc64:testcrc\tintact-crc:3E26AC3853066993\tfalse\tbinding site:21-21\tbinding site:24-24\t1\t1\tpsi-mi:\"MI:0078\"(nucleotide sequence identification)\tpsi-mi:\"MI:0078\"(nucleotide sequence identification)\n";
        Row row = mitabLineMapper.mapLine(mitab27, 0);

        Assert.assertNotNull(row);
        Assert.assertEquals(43, row.getAllFields().size());
    }

    @Test
    public void read_mitab_25_default_no_arguments() throws Exception {

        MitabCalimochoLineMapper mitabLineMapper = new MitabCalimochoLineMapper();

        String mitab25 = "uniprotkb:P73045\tintact:EBI-1579103\tintact:EBI-1607518\tintact:EBI-1607516\tuniprotkb:slr1767(locus name)\tuniprotkb:alias2(gene name)\tpsi-mi:\"MI:0018\"(two hybrid)\tauthor et al.(2007)\tpubmed:18000013\ttaxid:4932(yeasx)|taxid:4932(\"Saccharomyces cerevisiae (Baker's yeast)\")\ttaxid:1142(9sync)|taxid:1142(Synechocystis)\tpsi-mi:\"MI:0915\"(physical association)\tpsi-mi:\"MI:0469\"(IntAct)\tintact:EBI-1607514\tauthor-score:C\n";
        Row row = mitabLineMapper.mapLine(mitab25, 0);

        Assert.assertNotNull(row);
        // 15 MITAB fields + two more taxid fields
        Assert.assertEquals(17, row.getAllFields().size());
    }

    @Test
    public void read_mitab_26() throws Exception {

        MitabCalimochoLineMapper mitabLineMapper = new MitabCalimochoLineMapper();

        String mitab26 = "uniprotkb:P73045\tintact:EBI-1579103\tintact:EBI-1607518\tintact:EBI-1607516\tuniprotkb:slr1767(locus name)\tuniprotkb:alias2(gene name)\tpsi-mi:\"MI:0018\"(two hybrid)\tauthor et al.(2007)\tpubmed:18000013\ttaxid:4932(yeasx)|taxid:4932(\"Saccharomyces cerevisiae (Baker's yeast)\")\ttaxid:1142(9sync)|taxid:1142(Synechocystis)\tpsi-mi:\"MI:0915\"(physical association)\tpsi-mi:\"MI:0469\"(IntAct)\tintact:EBI-1607514\tauthor-score:C\t-\tpsi-mi:\"MI:0499\"(unspecified role)\tpsi-mi:\"MI:0499\"(unspecified role)\tpsi-mi:\"MI:0496\"(bait)\tpsi-mi:\"MI:0498\"(prey)\tpsi-mi:\"MI:0326\"(protein)\tpsi-mi:\"MI:0326\"(protein)\trefseq:NP_440386.1\trefseq:test2\tirefindex:test3\tcomment:test1\tcomment:test2\tcomment:test3\ttaxid:9606(human)\tkd:1.5\t2008/01/14\t2008/09/22\tcrc64:9E0E98F314F90177\tcrc64:testcrc\tintact-crc:3E26AC3853066993\tfalse\n";
        Row row = mitabLineMapper.mapLine(mitab26, 0);

        Assert.assertNotNull(row);
        // 36 fields -1 (no complex expansion) + 2 taxids more
        Assert.assertEquals(37, row.getAllFields().size());

    }
}
