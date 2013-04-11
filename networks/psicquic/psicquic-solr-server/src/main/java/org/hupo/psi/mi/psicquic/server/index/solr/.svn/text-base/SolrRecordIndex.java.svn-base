package org.hupo.psi.mi.psicquic.server.index.solr;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # SolrRecordIndex: access to solr-based index
 #         NOTE: supports only one set of solr server connections (base or 
 #               shards) shared between document submitting instances
 #   
 # ========================================================================== */

import org.apache.solr.common.SolrDocumentList;
import org.apache.solr.common.params.*;

import java.net.MalformedURLException;

import java.io.*;
import java.util.*;
import java.util.zip.*;

import java.util.regex.*;

import org.apache.solr.core.*;
import org.apache.solr.common.*;
import org.apache.solr.client.solrj.*;
import org.apache.solr.client.solrj.impl.*;
import org.apache.solr.client.solrj.response.*;

import org.hupo.psi.mi.psicquic.util.JsonContext;

import org.hupo.psi.mi.psicquic.server.index.RecordIndex;
import org.hupo.psi.mi.psicquic.server.*;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.lang.Thread;

public class SolrRecordIndex implements RecordIndex{
    
    PsqContext  psqContext = null;

    int queueSize = 4096;
    int solrconTCount = 2;

    private static SolrServer baseSolr = null;

    private String baseUrl = null;
    private String host = null;
    
    private String shardStr = "";
    
    private SolrParams defSolrParams = null;

    List<String> shardUrl = null;

    private static List<SolrServer> shSolr = null;

    CRC32 shCRC = null;

    private Map<String,Map> ricon = null;
    private Map<String,PsqTransformer> inTransformer = null;
    
    // record cache
    //-------------

    Map recordCache = null; 

    public Map getRecordCache(){
        return recordCache;
    }
    
    //--------------------------------------------------------------------------
    
    public void setPsqContext( PsqContext context ){
        this.psqContext = context;
        initialize( true );
    }

    public PsqContext getPsqContext(){
        return psqContext;
    }
    
    //--------------------------------------------------------------------------

    public SolrRecordIndex(){}
    
    public SolrRecordIndex( PsqContext context, int solrconTCount ){
        this.psqContext = context;
        this.solrconTCount = solrconTCount;
        initialize( true );
    }
    
    public SolrRecordIndex( PsqContext context, String host,
                            int solrconTCount ){
        this.psqContext = context;
        this.solrconTCount = solrconTCount;

        if( host != null ){
            this.host = host;
        }
        initialize( true );
    }
    
    public SolrRecordIndex( PsqContext context ){
        this.psqContext = context;
        initialize( true );
    }
    
    public SolrRecordIndex( PsqContext context, String host ){
        this.psqContext = context;

        if( host != null ){
            this.host = host;
        }
        initialize( true );
    }

    //--------------------------------------------------------------------------    

    public void initialize(){
        initialize( false );
    }
    
    public void initialize( boolean force ){
       
        // initialize cache
        //-----------------

        if( recordCache == null ){

            WeakHashMap whm = 
                new WeakHashMap( new HashMap<String,SolrDocument>());
            recordCache = Collections.synchronizedMap( whm ); 
        }

        if( force || baseUrl == null ){
            
            Log log = LogFactory.getLog( this.getClass() );
            log.info( " initilizing SolrRecordIndex: psqContext=" + psqContext );
	    if( psqContext != null ){
		log.info( "                        "
			  + "JsonConfig=" + psqContext.getJsonConfig() );
            }
            
            if( psqContext != null && psqContext.getJsonConfig() != null ){

                // SOLR URLs
                //----------

                Map jSolrCon = (Map) 
                    ((Map) psqContext.getJsonConfig().get("index")).get("solr");
                
                log.info( "    url=" + jSolrCon.get("url") );
                log.info( "    shard=" + jSolrCon.get("shard") );
                
                baseUrl = (String) jSolrCon.get("url");
                
                if( host != null ){
                    baseUrl = hostReset( baseUrl , host );
                }
                
                if( jSolrCon.get("shard") != null ){
                    List shardList = (List) jSolrCon.get("shard");
                    //String shardStr= "";
                    for( Iterator is = shardList.iterator(); is.hasNext(); ){
                        String csh = (String) is.next();
                        shardStr += csh + ",";
                        if( shardUrl == null ){
                            shardUrl = new ArrayList<String>();
                        }
                        //shardUrl.add( "http://" + csh );

                        if( host != null ){
                            shardUrl.add( hostReset( "http://" + csh, host ) );
                        } else {
                            shardUrl.add( "http://" + csh );
                        }
                    }
                    if( ! shardStr.equals("") ){
                        shardStr = shardStr.substring( 0, shardStr.length()-1);
                    }
                }
                
                // query parameters
                //-----------------
                
                defSolrParams 
                    = initSolrParams( (Map) jSolrCon.get( "param-default" ) );
                
                log.info( " initilizing SolrRecordIndex:"
                          + " param-default=" + defSolrParams );

                // transformations
                //----------------
                
                Map jFmtCon = (Map) ((Map) psqContext.getJsonConfig()
                                     .get("index")).get("transform");
                
                ricon = new HashMap<String,Map>();
                inTransformer = new HashMap<String,PsqTransformer>();

                for( Iterator ik = jFmtCon.keySet().iterator(); 
                     ik.hasNext(); ){
                    
                    String rtp = (String) ik.next();
                    Map rtr = (Map) jFmtCon.get( rtp );
                    
                    ricon.put( rtp, rtr );
                }
            }
        }
    }
    
    //--------------------------------------------------------------------------

    public String getBaseUrl(){
        return baseUrl;
    }
    
    //--------------------------------------------------------------------------

    public String getShardStr(){
        return shardStr;
    }
    
    //--------------------------------------------------------------------------

    private String hostReset( String url, String newHost ){

        if( host == null ) return url;

        // http://aaa:8888/d/d/d

        try{
            Pattern p = Pattern.compile( "([^/]*//)([^/:]+)(:[0-9]+)?(.*)" );
            Matcher m = p.matcher( url );

            if( m.matches() ){

                String prefix = m.group( 1 );
                String host = m.group( 2 );
                String port = m.group( 3 );
                String postfix = m.group( 4 );

                String newUrl = prefix + newHost + port + postfix;

                Log log = LogFactory.getLog( this.getClass() );
                log.info( "hostReset: old=" + url + " new=" + newUrl );

                return newUrl;

            } else {
                return url;
            }
        } catch(Exception ex ){
            return url;
        }
    }
    
    //--------------------------------------------------------------------------

    public ResultSet query( String query ){
        return query( query, null, -1, -1 );
    }

    public ResultSet query( String query, long firstRecord, long blockSize ){
        return query( query, null, firstRecord, blockSize );
    }
    
    public ResultSet query( String query, Map<String,List<String>> xquery,
                            long firstRecord, long blockSize ){
        
	Log log = LogFactory.getLog( this.getClass() );
        log.debug( "query: DefSolrParams=" + defSolrParams );
        
        ModifiableSolrParams params = new ModifiableSolrParams( defSolrParams );
        
        if( baseUrl == null ){ initialize(); }
        log.debug( "   SolrRecordIndex: baseUrl=" + baseUrl );

        if( baseUrl != null ){
                
            SolrServer solr = new HttpSolrServer( baseUrl );                       
                
            // set shards when needed
            //-----------------------
            
            if( shardStr != null ){
                params.set( "shards", shardStr );
            }
            
            // Query
            //------
            
            params.set( "q", query );
            
            // Paging
            //-------

            if( firstRecord >= 0 ){
                params.set( "start", Long.toString( firstRecord ) );
            }

            if( blockSize >= 0 ){
                params.set( "rows", Long.toString( blockSize ) );
            }

            // Extended query paramaters
            //--------------------------
            
            // -faceting: MiqlxGroupBy
            //------------------------
            
            if( xquery != null  && xquery.get( "MiqlxGroupBy:") != null ){
                List<String> ff = xquery.get( "MiqlxGroupBy:" );
                
                params.set( "facet", "true" );
                for( Iterator<String> fi = ff.iterator(); fi.hasNext(); ){
                        
                    // NOTE: restrict fields ???
                    
                    params.set( "facet.field", fi.next() + "_s" );
                }
            }
            
            // -faceting: FacetField (same as GroupBy ) 
            //-----------------------------------------

            if( xquery != null  && xquery.get( "MiqlxFacetField:") != null ){
                List<String> ff = xquery.get( "MiqlxFacetField:" );
                
                params.set( "facet", "true" );
                for( Iterator<String> fi = ff.iterator(); fi.hasNext(); ){
                        
                    // NOTE: restrict fields ???
                    // NOTE: 'set' limits faceting to one field only

                    params.set( "facet.field", fi.next() + "_s" );
                }
            }

            // -filter queries
            //----------------

            if( xquery != null  && xquery.get( "MiqlxFilter:") != null ){
                List<String> fq = xquery.get( "MiqlxFilter:" );

                for( Iterator<String> fi = fq.iterator(); fi.hasNext(); ){

                    // NOTE: multiple filters allowed

                    if( params.get( "fq" ) == null ){
                        params.set( "fq", fi.next() );
                    } else {
                        params.add( "fq", fi.next() );
                    }
                }
            }
            
            try{            
                QueryResponse response = solr.query( params );
                log.debug( "\n\nresponse:\n\n\n" + response +"\n\n\n");
                
                SolrDocumentList sdl = response.getResults();
                
                ResultSet rs = new ResultSet( sdl.getStart(),
                                              sdl.getNumFound(),
                                              (List) sdl );
                
                String recid = psqContext.getRecId();
                
                for( Iterator<SolrDocument> i = sdl.iterator(); i.hasNext(); ){
                    
                    SolrDocument cdoc = i.next();
                    String rec = (String) cdoc.get( recid );
                    
                    synchronized( recordCache ){
                        recordCache.put( rec, cdoc );
                    }
                }

                // Extended query results
                //-----------------------

                // - faceting
                //-----------

                if( response.getFacetFields() != null 
                    && response.getFacetFields().size() > 0 ){

                    log.info("processing facets");
                
                    Map<String,List<ValueCount>> rsFacets
                        = new HashMap<String,List<ValueCount>>();

                    for( Iterator<FacetField> iff = response.getFacetFields()
                             .iterator(); iff.hasNext(); ){

                        FacetField ff = iff.next();
                        String ffName = ff.getName().replaceAll( "_s$",""); 
                   
                        if( rsFacets.get( ffName )  == null ){
                            rsFacets.put( ffName, new ArrayList<ValueCount>());
                        }
     
                        for( Iterator<FacetField.Count> ifv = ff.getValues()
                                 .iterator(); ifv.hasNext(); ){
                            
                            FacetField.Count fv = ifv.next();
                            rsFacets.get( ffName )
                                .add( new ValueCount( fv.getName(), fv.getCount() ) );
                        }                        
                    }
                    
                    log.debug( "processing facets: done" );
                    rs.setMeta( "groups", rsFacets );
                    log.debug( "processing facets: meta=" + rs.getMeta() );
                }

                return rs;
                
            } catch( SolrServerException sex ){
                sex.printStackTrace();
            }
        }
        return null;
    }

    //--------------------------------------------------------------------------
    
    public void connect() throws MalformedURLException{

        Log log = LogFactory.getLog( this.getClass() );
        log.info( "SolrRecordIndex: connecting" );
        
        if( shardUrl.size() == 0 && baseSolr == null ){
	    
            baseSolr = new ConcurrentUpdateSolrServer( baseUrl, 
                                                       queueSize, 
                                                       solrconTCount );	    
        } else {
            log.info( "SolrRecordIndex: baseSolr already initialized: " 
                      + baseSolr );
        }
        
        if( shardUrl.size() > 0 && shSolr == null ){
            
            shSolr = new ArrayList<SolrServer>();
            
            for( Iterator<String> is = shardUrl.iterator(); is.hasNext(); ){
                
                shSolr.add( new ConcurrentUpdateSolrServer( is.next(), 
                                                            queueSize, 
                                                            solrconTCount ) );                  
	    }
        } else {
            if( shardUrl.size() > 0 ){
                log.info( "SolrRecordIndex: Shards already initialized: " 
                          + shSolr );
            }                           
        }
        
        log.info( "SolrRecordIndex: connected" );
    }

    //--------------------------------------------------------------------------

    public void add( SolrInputDocument doc ) 
        throws SolrServerException,IOException{
        
        baseSolr.add( doc );
    }
    
    //--------------------------------------------------------------------------
    
    public void add( int shard, SolrInputDocument doc ) 
        throws SolrServerException,IOException{
        
        shSolr.get( shard ).add( doc );
    }

    //--------------------------------------------------------------------------
    
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
    
    public void clear(){
        
        Log log = LogFactory.getLog( this.getClass() );
        
        try{
            if( shSolr != null && shSolr.size()>0 ){
                int i =0;
                for( Iterator<SolrServer> is = shSolr.iterator(); 
                     is.hasNext(); ){
                    
                    SolrServer css = is.next();
                    css.deleteByQuery( "*:*" );
                    css.commit();

                    log.info( "SolrRecordIndex Shard#"+i+" Cleared" );
                    i++;
                } 
            } else {
                baseSolr.deleteByQuery( "*:*" );
                baseSolr.commit();
                log.info( "SolrRecordIndex Base Index Cleared");
            }
        }catch( SolrServerException ssx ){
            log.error( "SolrRecordIndex Clear Error:", ssx );
        }catch( IOException  iox ){
            log.error( "SolrRecordIndex Clear Error:", iox );
        }
    }
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    public void addFile( File f, String name, String format, String compress){
        
	Log log = LogFactory.getLog( this.getClass() );
        Map rtr = ricon.get( format );
        
        List<SolrInputDocument> output = new ArrayList<SolrInputDocument>();

        // Known record transformers (add more when implemented)
        //------------------------------------------------------
        
        if( inTransformer.get( format ) == null ){
            if( ((String) rtr.get("type")).equalsIgnoreCase( "XSLT" ) ){
        
                log.info( " Initializing transformer: format=" + format
                          + " type=XSLT config=" + rtr.get( "config" ) );
                
                PsqTransformer recTr = 
                    new XsltTransformer( (Map) rtr.get("config") );
                
                inTransformer.put( format, recTr );
                log.info( " Initializing transformer(" + format
                          + "): new=" + inTransformer.get( format ) );
            }    
        
            if( ((String) rtr.get("type")).equalsIgnoreCase( "CALIMOCHO" ) ){
            		
                log.info( " Initializing transformer: format=" + format
                          + " type=CALIMOCHO config=" + rtr.get("config") );
                
                PsqTransformer recTr = 
                    new CalimochoTransformer( (Map) rtr.get( "config" ) );
                inTransformer.put( format, recTr );
                log.info( " Initializing transformer(" + format 
                          + "): new=" + inTransformer.get( format ) );                
            }
        }    

        PsqTransformer rt = (PsqTransformer) inTransformer.get( format );
        
        if( rt == null ){
            // NOTE: throw unknown transformer exception ?
        }
        
        try{
            if( compress!= null
                && compress.equalsIgnoreCase("zip") ){
                
                processZipFile( rt, name, new ZipFile( f ));
            } else {
                processFile( rt, name, new FileInputStream( f ));
            }

        } catch( Exception ex ){
            log.info( ex.getMessage(), ex );
            return;
        }
    }

    //--------------------------------------------------------------------------

    private void processZipFile( PsqTransformer rt,
                                 String fileName,  ZipFile zf  )
        throws java.io.IOException{
        
        for( Enumeration zfe = zf.entries(); zfe.hasMoreElements(); ){
            
            ZipEntry ze = (ZipEntry) zfe.nextElement();
            
            if( !ze.isDirectory() ){
                InputStream is = zf.getInputStream( ze );
                processFile( rt,
                             fileName + "::" + ze.getName() , is );
            }
        }
    }
    
    //--------------------------------------------------------------------------
    
    private void processFile( PsqTransformer rt,
                              String fileName,  InputStream is ){

        Log log = LogFactory.getLog( this.getClass() );
        rt.start( fileName, is );
        
        log.info( " SolrRecordIndex(add): transformation start (" 
                  + fileName + ")..." + rt );

        while( rt.hasNext() ){
            
            Map cdoc= rt.next();
		
            String recId = (String) cdoc.get("recId");;
            SolrInputDocument doc 
                = (SolrInputDocument) cdoc.get("solr");
            log.debug( " SolrRecordIndex(add): recId=" + recId );
            log.debug( " SolrRecordIndex(add): SIDoc=" + doc );
            
            try{
                if( shSolr.size() > 1 ){
                    
                    int shard = this.shardSelect( recId );
                    int shMax = shSolr.size() - 1;
                    
                    log.debug( " SolrRecordIndex(add): recId=" + recId +
                               " shard= " + shard + " (max= " + shMax +")" );
                    
                    log.debug( " SolrRecordIndex(add): document START");
                    log.debug( doc );
                    log.debug( " SolrRecordIndex(add): document END");
                    
                    this.add( shard, doc );
                } else {
                    this.add( doc );
                }
            } catch(Exception ex){
                ex.printStackTrace();
            }                    
        }

        try{
            // NOTE: commit every n (configurable) records ?
            this.commit();
            log.info( " SolrRecordIndex(add): commit DONE");
        }catch( Exception ex ){
            ex.printStackTrace();
        }
    }
    
    //--------------------------------------------------------------------------

    public void addFile( String format, String fileName, InputStream is ){
        
	Log log = LogFactory.getLog( this.getClass() );
        
        Map rtr = ricon.get( format );
        
        List<SolrInputDocument> output = new ArrayList<SolrInputDocument>();

        // Known record transformers (add more when implemented)
        //------------------------------------------------------
        
        if( inTransformer.get( format ) == null ){
            if( ((String) rtr.get("type")).equalsIgnoreCase( "XSLT" ) ){
        
                log.info( " Initializing transformer: format=" + format
                          + " type=XSLT config=" + rtr.get( "config" ) );
                
                PsqTransformer recTr = 
                    new XsltTransformer( (Map) rtr.get("config") );
                
                inTransformer.put( format, recTr );
                log.info( " Initializing transformer(" + format
                          + "): new=" + inTransformer.get( format ) );
            }    
        
            if( ((String) rtr.get("type")).equalsIgnoreCase( "CALIMOCHO" ) ){
            		
                log.info( " Initializing transformer: format=" + format
                          + " type=CALIMOCHO config=" + rtr.get("config") );
                
                PsqTransformer recTr = 
                    new CalimochoTransformer( (Map) rtr.get( "config" ) );
                inTransformer.put( format, recTr );
                log.info( " Initializing transformer(" + format 
                          + "): new=" + inTransformer.get( format ) );                
            }
        }    
                
        PsqTransformer recordTransformer 
            = (PsqTransformer) inTransformer.get( format );
        
        if( recordTransformer == null ){
            // NOTE: throw unknown transformer exception ?
        }
        
        recordTransformer.start( fileName, is );
        log.info( " SolrRecordIndex(add): transformation start (" 
                  + fileName + ")..." + recordTransformer);
        while( recordTransformer.hasNext() ){
		
            Map cdoc= recordTransformer.next();
		
            String recId = (String) cdoc.get("recId");;
            SolrInputDocument doc 
                = (SolrInputDocument) cdoc.get("solr");
            log.debug( " SolrRecordIndex(add): recId=" + recId );
            log.debug( " SolrRecordIndex(add): SIDoc=" + doc );
            
            try{
                if( shSolr.size() > 1 ){
                    
                    int shard = this.shardSelect( recId );
                    int shMax = shSolr.size() - 1;
                    
                    log.debug( " SolrRecordIndex(add): recId=" + recId +
                               " shard= " + shard + " (max= " + shMax +")" );

                    log.debug( " SolrRecordIndex(add): document START");
                    log.debug( doc );
                    log.debug( " SolrRecordIndex(add): document END");
                    
                    this.add( shard, doc );
                } else {
                    this.add( doc );
                }
            } catch(Exception ex){
                ex.printStackTrace();
            }                    
        }
        
        log.info( " SolrRecordIndex(add): transformation DONE");

        try{
            // NOTE: commit every n (configurable) records ? 
            this.commit();
            log.info( " SolrRecordIndex(add): commit DONE");
        }catch( Exception ex ){
            ex.printStackTrace();
        }
    }

    //--------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------
    
    private SolrParams initSolrParams( Map map ){
        
        if( map == null ) return new ModifiableSolrParams();
        
        Map<String,String[]> mmap = new HashMap<String,String[]>();
        
        for( Iterator<Map.Entry> mi = map.entrySet().iterator(); 
             mi.hasNext(); ){

            Map.Entry me = mi.next();

            String key = (String) me.getKey();
            Object val = me.getValue();
           
            try{  // NOTE: Some ugly casting; arrgh... ;o)

                if( val instanceof String ){

                    String[] va = new String[1]; 
                    va[0] = (String) val;
                    mmap.put( key, va );
                } else {
                    String[] vl = 
                        (String[]) ((List) val).toArray( new String[0] );
                    mmap.put( key, vl );
                }
            }catch( Exception ex ){
                Log log = LogFactory.getLog( this.getClass() );
                log.info( " initilizing SolrRecordIndex: " 
                          + "param-default format error=" + val );
                
                ex.printStackTrace();
            }
        }
        return new ModifiableSolrParams( mmap );
    }
}
