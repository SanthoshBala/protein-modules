package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id:: ValueCount.java 1130 2012-07-17 16:56:03Z lukasz99                    $
 # Version: $Rev:: 1130                                                        $
 #==============================================================================
 #                                                                             
 # Value.Count: utility class
 #                                                                             
 #=========================================================================== */

public class ValueCount{

    private String value;
    private long count;
    
    public ValueCount( String value, long count){
        this.value = value;
        this.count = count;
    }

    public String getValue(){
        return value;
    }

    public long getCount(){
        return count;
    }
}
