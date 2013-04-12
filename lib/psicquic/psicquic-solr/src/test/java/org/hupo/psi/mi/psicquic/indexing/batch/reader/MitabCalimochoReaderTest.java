package org.hupo.psi.mi.psicquic.indexing.batch.reader;

import junit.framework.Assert;
import org.hupo.psi.calimocho.model.Row;
import org.junit.Test;
import org.springframework.core.io.FileSystemResource;

/**
 * Unit tester of the MitabCalimocho reader
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>12/07/12</pre>
 */

public class MitabCalimochoReaderTest {

    @Test
    public void read_mitab_27() throws Exception {

        MitabCalimochoReader reader = new MitabCalimochoReader();
        reader.setResource(new FileSystemResource(MitabCalimochoReader.class.getResource("/samples/sampleFileNegative.txt").getFile()));

        reader.doOpen();

        Row row = reader.read();

        reader.close();

        Assert.assertNotNull(row);
        Assert.assertEquals(139, row.getAllFields().size());
    }

    @Test
    public void read_mitab_25() throws Exception {

        MitabCalimochoReader reader = new MitabCalimochoReader();
        reader.setResource(new FileSystemResource(MitabCalimochoReader.class.getResource("/samples/sampleFileMitab25.txt").getFile()));

        reader.doOpen();

        Row row = reader.read();

        reader.close();

        Assert.assertNotNull(row);
        Assert.assertEquals(28, row.getAllFields().size());
    }

    @Test
    public void read_mitab_26() throws Exception {

        MitabCalimochoReader reader = new MitabCalimochoReader();
        reader.setResource(new FileSystemResource(MitabCalimochoReader.class.getResource("/samples/sampleFileMitab26.txt").getFile()));

        reader.doOpen();

        Row row = reader.read();

        reader.close();

        Assert.assertNotNull(row);
        Assert.assertEquals(125, row.getAllFields().size());

    }
}
