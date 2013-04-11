package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #                                                                             
 # ResultSet: index query results
 #                                                                             
 #=========================================================================== */

import java.util.*;

public class ResultSet{

    List resultList;

    long firstResult;
    long maxResult;

    String format;

    Map<String,Object> meta;
    
    public ResultSet(){
        firstResult = 0;
        maxResult = 0;
    }

    public ResultSet( String format ){
        firstResult = 0;
        maxResult = 0;
        this.format = format;
    }

    public ResultSet( long first, long max, List results ){
        resultList = results;
        firstResult = first;
        maxResult = max;
    }

    public ResultSet( String format, long first, long max, 
                      List results ){
        resultList = results;
        firstResult = first;
        maxResult = max;
        this.format = format;
    }

    public List getResultList(){
        if(resultList == null ){
            resultList = new ArrayList();
        }
        return resultList;
    }

    public void setResultList( List results ){
        resultList = results;
    }
    
    public long getFirstResult(){
        return firstResult;
    }

    public void setFirstResult( long first ){
        firstResult = first;
    }

    public long getMaxResult(){
        return maxResult;
    }

    public void setMaxResult( long max ){
        maxResult = max;
    }

    public String getFormat(){
        return format;
    }
    
    public Map<String,Object> getMeta(){
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