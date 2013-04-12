package org.hupo.psi.mi.psicquic.server.index.solr;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # SolrIndex: access to solr-based index
 #
 #=========================================================================== */

import org.apache.solr.common.params.ModifiableSolrParams;
import org.apache.solr.common.SolrDocumentList;

import java.net.MalformedURLException;

import java.io.*;
import java.util.*;
import java.util.zip.CRC32;

import org.apache.solr.core.*;
import org.apache.solr.common.*;
import org.apache.solr.client.solrj.*;
import org.apache.solr.client.solrj.impl.*;
import org.apache.solr.client.solrj.response.*;

import org.hupo.psi.mi.psicquic.util.JsonContext;

import org.hupo.psi.mi.psicquic.server.index.*;
import org.hupo.psi.mi.psicquic.server.*;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class SolrIndex implements Index{

    JsonContext context = null;

    String baseUrl = null;
    SolrServer baseSolr = null;

    List<String> shardUrl = null;
    List<SolrServer> shSolr = null;
    CRC32 shCRC = null;

    Map<String,Map> ricon = null;
    Map<String,Map> rscon = null;
    
    public void setContext( JsonContext context ){
        this.context = context;
    }

    public SolrIndex(){}
    
    public SolrIndex( JsonContext context ){
        this.context = context;
        initialize( true );
    }
    
    public void initialize(){
        initialize( false );
    }
    
    public void initialize( boolean force ){
        
        if( force || baseUrl == null ){
            
            Log log = LogFactory.getLog( this.getClass() );
            log.info( " initilizing SolrIndex" );
            
            if( context != null && context.getJsonConfig() != null ){

                // SOLR URLs
                //----------

                Map jSolrCon = (Map) 
                    ((Map) context.getJsonConfig().get("index")).get("solr");

                log.info( "    url=" + jSolrCon.get("url") );
                log.info( "    shard=" + jSolrCon.get("shard") );
                
                baseUrl = (String) jSolrCon.get("url");
                
                if( jSolrCon.get("shard") != null ){
                    List shardList = (List) jSolrCon.get("shard");
                    String shardStr= "";
                    for( Iterator is = shardList.iterator(); is.hasNext(); ){
                        String csh = (String) is.next();
                        shardStr += csh + ",";
                        if( shardUrl == null ){
                            shardUrl = new ArrayList<String>();
                        }
                        shardUrl.add( "http://" + csh );
                    }
                    if( !shardStr.equals("") ){
                        shardStr = shardStr.substring( 0, shardStr.length()-1);
                    }
                }
                
                Map jFmtCon = (Map) ((Map) context.getJsonConfig().get("index"))
                    .get("transform");
                
                ricon = new HashMap<String,Map>();

                for( Iterator ik = jFmtCon.keySet().iterator(); ik.hasNext(); ){
                    
                    String rtp = (String) ik.next();
                    Map rtr = (Map) jFmtCon.get( rtp );

                    ricon.put( rtp, rtr );
                }
            }
        }
    }
    
    //--------------------------------------------------------------------------
    
    public ResultSet query( String query ){

        ResultSet rs = new ResultSet();
	
        try{
            if( baseUrl == null ){ initialize(); }
            if( baseUrl != null ){

                SolrServer solr = new CommonsHttpSolrServer( baseUrl );
                ModifiableSolrParams params = new ModifiableSolrParams();
                params.set( "q", query );
                
                // set shards if needed ?
                //-----------------------
                
                try{            
                    QueryResponse response = solr.query( params );
                    System.out.println( "response = " + response );

                    SolrDocumentList res = response.getResults();
                    rs.setResultList( res );
                    
                } catch( SolrServerException sex ){
                    sex.printStackTrace();
                }
            }
        } catch( MalformedURLException mex ){
            mex.printStackTrace();
        }
        return rs;
    }

    //--------------------------------------------------------------------------
    
    public void connect() throws MalformedURLException{

        
        baseSolr = new CommonsHttpSolrServer( baseUrl );
        
        Log log = LogFactory.getLog( this.getClass() );
        log.info( "   SolrIndex: connecting" );

        
        if( shardUrl.size() > 0 ){

            shSolr = new ArrayList<SolrServer>();
            
            for( Iterator<String> is = shardUrl.iterator(); is.hasNext(); ){
                shSolr.add( new CommonsHttpSolrServer( is.next() ) );                        
            }
        } 

        log.info( "              shards: total=" + shardUrl.size() 
                  + " connected=" + shSolr.size() );

    }
    
    public void add( SolrInputDocument doc ) 
        throws SolrServerException,IOException{
        
        baseSolr.add( doc );
    }
    
    public void add( int shard, SolrInputDocument doc ) 
        throws SolrServerException,IOException{
        
        shSolr.get( shard ).add( doc );
    }

    
    public void commit()
        throws SolrServerException,IOException{

        if( shSolr != null && shSolr.size() > 0 ){
            for( Iterator<SolrServer> is = shSolr.iterator(); is.hasNext(); ){
                is.next().commit();
            } 
        } else {
            baseSolr.commit();
        }    
    }

    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    public void addFile( String format, String fileName, InputStream is ){
        
        Map rtr = ricon.get( format );
        
        List<SolrInputDocument> output = new ArrayList<SolrInputDocument>();

        if( ((String) rtr.get("type")).equalsIgnoreCase("XSLT") ){
            if( rtr.get( "transformer" ) == null ){

                Log log = LogFactory.getLog( this.getClass() );
                log.info( " Initializing transformer: format=" + format
                          + " type=XSLT config=" + rtr.get("config") );

                PsqTransformer rt = 
                    new XsltTransformer( (String) rtr.get("config") );
                rtr.put( "transformer",rt );
            }
            
            PsqTransformer rt = (PsqTransformer) rtr.get( "transformer" );
            rt.start( fileName, is );
            
	    while( rt.hasNext() ){

		SolrInputDocument doc = rt.next();
		SolrInputField pidf = doc.getField("pid");
		String pid = (String) pidf.getFirstValue();
		
		try{
		    if( shSolr.size() > 1 ){
			
			int shard = this.shardSelect( pid );
			Log log = LogFactory.getLog( this.getClass() );
			int shMax = shSolr.size() - 1;
			log.info( " SolrIndex(add): pid=" + pid 
				  + " shard= " + shard 
				  + " (max= " + shMax +")" );
			
			this.add( shard, doc );
		    } else {
			this.add( doc );
		    }
		} catch(Exception ex){
		    ex.printStackTrace();
		}                    
	    }
            
            try{
                this.commit();   
            }catch( Exception ex ){
                ex.printStackTrace();
            }
        }   
    }

    public int shardSelect( String id ){
        
        if( shCRC == null ){
            shCRC = new CRC32();
        }

        shCRC.reset();
        shCRC.update( id.getBytes() ); 
        
        Long crc =  new Long( shCRC.getValue() );
        int shc = shSolr.size();
      
        if( crc.intValue() > 0 ){
            return crc.intValue() % shc;
        } else {
            return - crc.intValue() % shc;
        }
    }

}