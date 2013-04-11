package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
# $Id::                                                                       $
# Version: $Rev::                                                             $
#==============================================================================
#                                                                             
# PsqResultSet: 
#                                                                             
#=========================================================================== */

import java.util.*;

public class PsqRecordSet{

    List recordList;

    int firstRecord;
    int maxRecord;
    
    Map<String,Object> meta;
    
    public PsqRecordSet(){
        recordList = new ArrayList();
        firstRecord = 0;
        maxRecord = 0;
    }
    
    public PsqRecordSet( int first, int max, List records ){
        recordList = records;
        firstRecord = first;
        maxRecord = max;
    }

    public List getRecordList(){
        return recordList;
    }

    public void setRecordList( List records ){
        recordList = records;
    }
    
    public int getFirstRecord(){
        return firstRecord;
    }

    public void setFirstRecord(int first ){
        firstRecord = first;
    }

    public int getMaxRecord(){
        return maxRecord;
    }

    public void setMaxRecord( int max ){
        maxRecord = max;
    }
    
    public Map<String,Object> getMeta(){

        if( meta == null ){
            meta = new HashMap<String,Object>();
        }
        return meta;
    }

    public void setMeta( Map<String,Object> map ){
        meta = map;
    }

    public void setMeta( String key, Object value ){

        if( meta == null ){
            meta = new HashMap<String,Object>();
        }

        meta.put( key, value );
    }

}