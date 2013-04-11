package org.hupo.psi.mi.psicquic.indexing.batch.reader;

import org.hupo.psi.calimocho.model.Row;
import org.hupo.psi.calimocho.tab.io.DefaultRowReader;
import org.hupo.psi.calimocho.tab.io.RowReader;
import org.hupo.psi.calimocho.tab.model.ColumnBasedDocumentDefinition;
import org.hupo.psi.calimocho.tab.util.MitabDocumentDefinitionFactory;
import org.springframework.batch.item.file.LineMapper;

/**
 * Mapper that will map a MITAB line to a Calimocho Row
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>29/05/12</pre>
 */

public class MitabCalimochoLineMapper implements LineMapper<Row>{

    private ColumnBasedDocumentDefinition documentDefinition;
    private RowReader rowReader;

    public MitabCalimochoLineMapper() {

        // if no version, then we suppose it is MITAB 25
        documentDefinition = MitabDocumentDefinitionFactory.mitab27();

        this.rowReader = new DefaultRowReader( documentDefinition );
    }

    /**
     * Convert a MITAB line in a calimocho row
     * @param line : mitab line
     * @param lineNumber
     * @return the matching Calimocho Row
     * @throws Exception
     */
    public Row mapLine(String line, int lineNumber) throws Exception {
        try {

            if (line == null){
                return null;
            }

            Row row = rowReader.readLine(line);

            return row;
        } catch (Exception e) {
            throw new Exception("Problem converting to binary interaction line "+lineNumber+": "+line, e);
        }
    }
}
