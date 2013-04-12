package org.hupo.psi.mi.psicquic.indexing.batch.reader;

import org.hupo.psi.calimocho.model.Row;
import org.springframework.batch.item.file.FlatFileItemReader;

/**
 * Reader of MITAB file
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>29/05/12</pre>
 */

public class MitabCalimochoReader extends FlatFileItemReader<Row> {

    public MitabCalimochoReader(){
        setLineMapper(new MitabCalimochoLineMapper());
    }

    @Override
    protected void doOpen() throws Exception {

        super.doOpen();
    }
}
