package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # PsicquicServer: 
 #
 #=========================================================================== */

import java.util.*;

import java.io.InputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.hupo.psi.mi.psicquic.*;
import org.hupo.psi.mi.psicquic.server.index.*;
import org.hupo.psi.mi.psicquic.server.store.*;

import org.hupo.psi.mi.psicquic.util.JsonContext;

public class PsicquicServer {
    
    PsqContext psqContext;
    
    public void setPsqContext( PsqContext context ){
        psqContext = context;
    }

    public PsqContext getPsqContext(){
        return psqContext;
    }
    
    //--------------------------------------------------------------------------

    private void initialize() {
        initialize( false );
    }
    
    //--------------------------------------------------------------------------
    
    private void initialize( boolean force) {

	Log log = LogFactory.getLog( this.getClass() );
	log.info( " psqContext=" + psqContext );
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    public ResultSet getByQuery( String query, String resultType,
                                 long firstResult, long blockSize ){
        
        Log log = LogFactory.getLog( this.getClass() );
        log.info( "getByQuery: context =" + psqContext);
        log.info( "getByQuery: q=" + query );
        
        log.info( " FR(in)=" + firstResult + " BS(in)=" + blockSize );
        
        // PAGING
        //-------
        
        if( firstResult < 0 ){
            firstResult = psqContext.getDefaultFirstResult();
        } 
        
        if( blockSize < 0 ){
            blockSize = psqContext.getDefaultBlockSize();
        } 

        // MIQLX
        //------
        
        Map<String,List<String>> miqlx = null;
        
        if( query != null && query.indexOf( " Miqlx" ) > -1 ){
            MiqlxFilter mf = new MiqlxFilter( psqContext );
            query = mf.process( query );
            miqlx = mf.getMiqlx();
        }
        
        log.info(":" + query +":");

        String viewType = psqContext.getDefaultView();
        
        if( miqlx != null ){
            for( Iterator mi = miqlx.entrySet().iterator(); mi.hasNext(); ){
                Map.Entry me = (Map.Entry) mi.next();

                log.info( "MIQLX: field=\'" + me.getKey() 
                          + "\'  value=\'" + me.getValue() +"\'");
            }

            if( miqlx.get("MiqlxView:") != null ){
                viewType = ((List<String>) miqlx.get("MiqlxView:")).get(0);
                

            }
        }
        log.info( " FR(q)=" + firstResult + " BS(q)=" + blockSize );
        ResultSet qrs = psqContext.getActiveIndex()
            .query( query, miqlx, firstResult, blockSize );
        
        ResultSet prs = 
            new ResultSet( viewType, qrs.getFirstResult(),
                           qrs.getMaxResult(), new ArrayList() );
        
        prs.setMeta( qrs.getMeta() );

	log.info( "getByQuery: rs="+ qrs); 
        
        for( Iterator i = qrs.getResultList().iterator(); i.hasNext(); ){
            Map in = (Map) i.next();
            log.debug( "getByQuery: in="+ in);

            String recId = (String) in.get( psqContext.getRecId() );

            String drecord =  psqContext.getActiveStore()
                .getRecord( recId , viewType );
        
            log.debug( " SolrDoc: recId=" + recId + " :: "  + drecord );
            prs.getResultList().add( drecord );
        }
        
        return prs;
    }

    //--------------------------------------------------------------------------

    public String toString( ResultSet rset, boolean records, boolean meta ){
        
        String rstr = "";

        // Records
        //--------

        if( records ){
            rstr += psqContext.getActiveStore().toString( rset );
        }
        
        // Meta info
        //----------

        if( meta && rset.getMeta() != null ){
            if( rset.getMeta().get( "groups") != null ){
                
                Map<String,List<ValueCount>> groups = 
                    (Map<String,List<ValueCount>>) rset.getMeta()
                    .get( "groups");
                
                for( Iterator<Map.Entry<String,List<ValueCount>>> ig =
                         groups.entrySet().iterator(); ig.hasNext(); ){
                    
                    Map.Entry<String,List<ValueCount>> me = ig.next();
                    
                    String field = me.getKey();
                    
                    for( Iterator<ValueCount> ivc = me.getValue().iterator();
                         ivc.hasNext(); ){

                        ValueCount vc = ivc.next();
                        String value = vc.getValue();
                        long count = vc.getCount();
                        rstr += "#" + "miqlx:groupby" + "\t"
                            + field + "\t"
                            + value + "\t"
                            + count + "\n";
                    }
                }
            }
        }
        return rstr;
    }
    
    //--------------------------------------------------------------------------

    public List<String> getSupportedReturnTypes( String service) {
        
        return (List<String>) ((Map) ((Map) psqContext.getJsonConfig()
                                      .get( "service" )).get( service ))
            .get( "supported-return-type" );
    }
    
    //--------------------------------------------------------------------------
    
    public String getVersion( String service ){
        
        return (String) ((Map) ((Map) psqContext.getJsonConfig()
                                .get( "service" )).get( service ))
            .get( "version" );
    }

    //--------------------------------------------------------------------------
    
    public List<String> getSupportedDbAcs( String service ){
        
        return (List<String>) ((Map) ((Map) psqContext.getJsonConfig()
                                      .get( "service" )).get( service ))
            .get( "supported-db-ac" );
    }
    
    //--------------------------------------------------------------------------

    public String getProperty( String service, String property ){     
        
        return (String) ((Map) ((Map) ((Map) psqContext.getJsonConfig()
                                       .get( "service" )) .get( service )) 
                         .get( "properties" ))
            .get( property );
    }
    
    //--------------------------------------------------------------------------
    
    public Set<Map.Entry> getProperties( String service ){
        
        Map propmap = (Map) ((Map) ((Map) psqContext.getJsonConfig()
                                    .get( "service" )).get( service ))
            .get( "properties" );
        
        return (Set<Map.Entry>) propmap.entrySet(); 
    }
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    public String buildQuery( String field, String db, String ac ){

        StringBuffer sb = new StringBuffer();
        sb.append( field ).append(":");
        if( db != null && !db.equals("") ){
            sb.append( db ).append( "\\:" );
        }
        sb.append( ac );
    
        return sb.toString();
    }

    //--------------------------------------------------------------------------
    
    public String buildQuery( String field, List<String> dbl, 
                              List<String> acl, String operand ){
        
        StringBuffer sb = new StringBuffer();
        sb.append( field ).append( ":(" );

        int pos = 0;
        for ( Iterator<String> idbl = dbl.iterator(); idbl.hasNext(); ){
            String db  = idbl.next();
            String ac = acl.get( pos );
            pos++;
            
            if( db != null && !db.equals("") ){
                sb.append( db ).append( "\\:" );
            }
            sb.append( ac );
            
            if( idbl.hasNext() ){
                sb.append( " " ).append( operand ).append( " " );
            }
        }

        sb.append( ")" );
        return sb.toString();
    }

    //--------------------------------------------------------------------------

    public String getHeader( String viewType ){
        return psqContext.getHeader( viewType );         
    }

    //--------------------------------------------------------------------------

    public String getFooter( String viewType ){
        return psqContext.getFooter( viewType );
    }
}
