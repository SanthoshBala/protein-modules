package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # MiqlxFilter: process MIQLX fields 
 #
 #=========================================================================== */

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class MiqlxFilter{
    
    PsqContext psqContext = null;

    public MiqlxFilter( PsqContext context ){
        psqContext = context;
    }
    
    Map<String, List<String>> miqlx = null;
    
    public Map<String, List<String>> getMiqlx(){
        return miqlx;
    }
    
    public String process( String query ){
        
        if( query == null || query.indexOf( " Miqlx" ) == -1 ) return query;

        List<String> miqlxFld = new ArrayList();
        
        List flst = (List) psqContext.getMiqlxDef().get( "field" );
        for( Iterator fi = flst.iterator(); fi.hasNext(); ){
            String  fld = (String) ((Map) fi.next()).get( "name" );
            miqlxFld.add( fld +":");
        }
        
        for( Iterator<String> mfi = miqlxFld.iterator(); mfi.hasNext(); ){
            
            if( query.indexOf( " Miqlx" ) == -1 ) return query;
            query = testMiqlx( query, mfi.next() );
        }        
        
        return query;
    }

    private String testMiqlx( String query, String field ){
        
        String nquery = query;

        Log log = LogFactory.getLog( this.getClass() );
        log.debug( " query=" + "<  field=" + field + "<" );
               

        while( query.indexOf( " " + field ) > 0 ){

            int vtStart = query.indexOf( " " + field );
            int vtStop = query.indexOf( ' ', vtStart + field.length() );
    
            String fval = null;
            
            if( vtStop > vtStart ){
                fval = query.substring( vtStart + field.length()+1, vtStop );
            } else {
                fval = query.substring( vtStart + field.length() + 1 );
            }
            
            nquery = query.substring( 0, vtStart );

            if( vtStop > vtStart ){
                nquery = nquery + query.substring( vtStop );
            }
            
            if( fval.length() > 0 ){
                if( miqlx == null ){
                    miqlx = new HashMap<String,List<String>>();
                }
                if( miqlx.get(field) == null ){
                    miqlx.put( field, new ArrayList<String>() );
                }
                miqlx.get( field ).add( fval );
            }
            query = nquery;
        }

        log.debug( " nquery=" + nquery + "<" );        
        return nquery;
    }
}
