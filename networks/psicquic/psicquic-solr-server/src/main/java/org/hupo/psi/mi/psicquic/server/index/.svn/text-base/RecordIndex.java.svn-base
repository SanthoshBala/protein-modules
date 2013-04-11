package org.hupo.psi.mi.psicquic.server.index;

import java.net.MalformedURLException;
import java.io.File;

import java.util.Map;
import java.util.List;

import org.hupo.psi.mi.psicquic.server.ResultSet;

public interface RecordIndex{

    public ResultSet query( String query );
    public ResultSet query( String query, long firstResult, long blockSize );
    public ResultSet query( String query, Map<String,List<String>> xquery,
                            long firstResult, long blockSize );
    public void initialize();
    public void clear();
    public void connect() throws MalformedURLException;
    //public void addFile( String format, String fileName, InputStream is );
    public void addFile( File f, String name, String format, String compress );
}
