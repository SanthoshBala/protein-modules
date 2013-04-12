package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #                                                                           
 # PsqTransformer: transforms input file into one or more solr input documents
 #                                                                              
 #=========================================================================== */

import org.apache.solr.common.SolrInputDocument;
import java.io.InputStream;

public interface PsqTransformer{

    public void start( String fileName, InputStream is );
    public void start( String fileName );

    public boolean hasNext();
    public SolrInputDocument next();

}
