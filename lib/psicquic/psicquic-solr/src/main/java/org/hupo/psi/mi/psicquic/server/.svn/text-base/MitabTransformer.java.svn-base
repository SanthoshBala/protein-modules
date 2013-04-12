package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # MitabTransformer: transform the incoming xml file into solr documents
 #
 #=========================================================================== */

import java.util.*;
import java.io.*;
import org.apache.solr.common.SolrInputDocument;

import org.hupo.psi.calimocho.model.Row;
import org.hupo.psi.calimocho.tab.io.DefaultRowReader;
import org.hupo.psi.calimocho.tab.model.ColumnBasedDocumentDefinition;
import org.hupo.psi.calimocho.tab.util.MitabDocumentDefinitionFactory;
import psidev.psi.mi.calimocho.solr.converter.Converter;

public class MitabTransformer implements PsqTransformer{

    private List<? extends Row> rows;
    int currIndex;
    
    public MitabTransformer( List<? extends Row> items ){
        this.rows = items;
        currIndex = -1;
    }

    public void start( String fileName ){
        try {
            start (fileName, new FileInputStream(new File(fileName)));
        } catch (Exception ex) {
            //TODO - do whatever you want to do with this exception
        }    
    }

    public void start( String fileName, InputStream is ){
        try {
            //TODO - check for mitab version
            ColumnBasedDocumentDefinition docDef = MitabDocumentDefinitionFactory.mitab25(); 
            DefaultRowReader reader = new DefaultRowReader(docDef);
            rows = reader.read( is );
        } catch (Exception ex) {
            //TODO - do whatever you want to do with this exception
        }
    }	

    public boolean hasNext(){
	return ( currIndex < rows.size() );
    }

    public SolrInputDocument next(){
        SolrInputDocument doc = null;
        currIndex++;
        try {
            Converter solrDocumentConverter = new Converter();
            doc =  solrDocumentConverter.toSolrDocument(rows.get(currIndex));
        } catch (Exception ex) {
            //do whatever you want to do with this exception
        }
        return doc;
    }
   
 
}