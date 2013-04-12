package org.hupo.psi.mi.psicquic.indexing.batch.writer;

import org.apache.solr.common.SolrInputDocument;
import org.hupo.psi.calimocho.model.Row;
import psidev.psi.mi.calimocho.solr.converter.Converter;

import java.util.Iterator;
import java.util.List;

/**
 * TODO comment this
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>23/08/12</pre>
 */

public class SolrDocumentCalimochoRowIterator implements Iterator<SolrInputDocument> {

    private Iterator<? extends Row> rowIterator;
    private Converter solrConverter;

    public SolrDocumentCalimochoRowIterator(List<? extends Row> items, Converter solrConverter){
        rowIterator = items.iterator();
        this.solrConverter = solrConverter;
    }

    public boolean hasNext() {
        return rowIterator.hasNext();
    }

    public SolrInputDocument next() {

        Row binaryInteraction = rowIterator.next();
        try {
            SolrInputDocument solrInputDoc = solrConverter.toSolrDocument(binaryInteraction);

            return solrInputDoc;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove");
    }
}
