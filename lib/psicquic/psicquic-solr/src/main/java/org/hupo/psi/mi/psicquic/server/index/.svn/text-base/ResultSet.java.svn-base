package org.hupo.psi.mi.psicquic.server.index;

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

    int firstResult;
    int maxResult;

    Map<String,Object> meta;
    
    public ResultSet(){
        resultList = new ArrayList();
        firstResult = 0;
        maxResult = 0;

    }
    
    public ResultSet( int first, int max, List results ){
        resultList = results;
        firstResult = first;
        maxResult = max;
    }

    public List getResultList(){
        return resultList;
    }

    public void setResultList( List results ){
        resultList = results;
    }
    
    public int getFirstResult(){
        return firstResult;
    }

    public void setFirstResult(int first ){
        firstResult = first;
    }

    public int getMaxResult(){
        return maxResult;
    }

    public void setMaxResult( int max ){
        maxResult = max;
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