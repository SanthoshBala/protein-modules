package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id:: PsqTransformer.java 1454 2012-09-26 22:54:34Z lukasz99                $
 # Version: $Rev:: 1454                                                        $
 #==============================================================================
 #                                                                           
 # PsqTransformer: transforms input file into one or more solr input documents
 #                                                                              
 #=========================================================================== */

import java.io.InputStream;
import java.util.Map;

public interface PsqTransformer{

    public void start( String fileName, InputStream is );
    public void start( String fileName );

    public boolean hasNext();
    public Map next();

    //public Map toString( RecordSet rset );

}
