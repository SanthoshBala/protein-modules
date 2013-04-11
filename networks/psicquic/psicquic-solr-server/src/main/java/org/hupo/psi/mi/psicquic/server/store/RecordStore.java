package org.hupo.psi.mi.psicquic.server.store;

/* =============================================================================
 # $Id:: RecordStore.java 1455 2012-09-28 01:15:16Z lukasz99                   $
 # Version: $Rev:: 1455                                                        $
 #==============================================================================
 #
 # RecordStore: record access 
 #
 #=========================================================================== */

import java.util.List;
import java.io.File;

import org.hupo.psi.mi.psicquic.server.*;

public interface RecordStore{
    public void initialize();
    public void clear();
    public void clearLocal();
    public void shutdown();
    public void addRecord( String id, String record, String format );
    //public void addFile( String fileName, String format, InputStream is );
    public void addFile( File f, String name, String format, String compress );

    public String getRecord( String id, String format );
    public List<String> getRecordList( List<String> id, String format );
    public String toString( ResultSet rset );
}
