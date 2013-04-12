package org.hupo.psi.mi.psicquic.server.builder;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # IndexBuilder: populate  PSICQUIC index 
 #
 #=========================================================================== */

import java.util.*;
import java.util.zip.*;
import java.net.*;
import java.io.*;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.hupo.psi.mi.psicquic.server.*;

import org.hupo.psi.mi.psicquic.server.index.*;
import org.hupo.psi.mi.psicquic.server.index.solr.SolrRecordIndex;

import org.hupo.psi.mi.psicquic.server.store.*;
import org.hupo.psi.mi.psicquic.server.store.derby.DerbyRecordStore;
import org.hupo.psi.mi.psicquic.server.store.solr.SolrRecordStore;

import org.hupo.psi.mi.psicquic.util.JsonContext;

import java.lang.Thread.State;

public class IndexBuilder{

    private List<String> solrURL = 
        Arrays.asList( "http://127.0.0.1:8080/psicquic-solr/ws/psq" );    

    private String rmgrURL = 
        "http://127.0.0.1:8080/psicquic-solr/recordmgr";
    
    public static final String 
        CONTEXT = "/etc/psq-context-default.json";
    
    public static final String 
        RFRMT = "mif254";

    Log log = null;
        
    String format = null;
    
    String root = null;
    String source = null;
    boolean zip = false;

    RecordIndex recordIndex = null;
    RecordStore recordStore = null;

    int builderTCount = 1;
    int solrconTCount = 1;
    

    String host = null;
    
    //--------------------------------------------------------------------------
    
    PsqContext psqContext = null;
        
    public PsqContext getPsqContext(){
        return psqContext;
    }
    
    //--------------------------------------------------------------------------
    
    public IndexBuilder( String ctx, String host, int btc, int stc, 
                         String format, boolean zip ){
        
        InputStream isCtx = null;
        this.log = LogFactory.getLog( this.getClass() );
        try{
            isCtx = new FileInputStream( ctx );
        } catch( Exception ex ){
            log.info( ex.getMessage(), ex );
        }
        
        this._IndexBuilder( isCtx, host, btc, stc, format, zip );          
    }
    
    public IndexBuilder( InputStream isCtx, String host, int btc, int stc,
                         String format, boolean zip ){

        this.log = LogFactory.getLog( this.getClass() );
        this._IndexBuilder( isCtx, host, btc, stc, format, zip );
    }

    private void _IndexBuilder( InputStream isCtx,  String host, 
                                int btc, int stc, String format, boolean zip ){
        this.zip = zip;
        this.format = format;
        this.source = source;
        
        log.info( " initilizing IndexBuilder: threads=" + builderTCount );
        
        // get context
        //------------

        JsonContext jsq = new JsonContext();
        
        try{                                
            jsq.readJsonConfigDef( isCtx );
        } catch( Exception ex ){
            log.info( ex.getMessage(), ex );
        }        

        psqContext = new PsqContext( jsq );
        
        Map jctx = psqContext.getJsonConfig();
        
        if( btc > 0 ){
            builderTCount = btc;
        } else {
            Integer btcInt = (Integer) ((Map)jctx.get( "builder" ))
                .get( "builder-thread-count" );           
            builderTCount = btcInt.intValue();
        }
        
        if( stc > 0 ){
            solrconTCount = stc;
        } else {
            Integer stcInt = (Integer) ((Map)jctx.get( "builder" ))
                .get( "solr-thread-count" );
            solrconTCount = stcInt.intValue();
        }

        if( host != null ){
            this.host = host;
        }
    }

    //--------------------------------------------------------------------------

    public void  clear(){

        Map jctx = (Map) psqContext.getJsonConfig();
         
        try{
            
            // clear index
            //------------

            String activeIndexName = psqContext.getActiveIndexName();
            
            if( activeIndexName.equals( "solr" ) ){
                recordIndex = new SolrRecordIndex( psqContext );
                recordIndex.initialize();
                recordIndex.connect();
                recordIndex.clear();
            }
            
            // clear data store
            //-----------------

            String activeStoreName = psqContext.getActiveStoreName();
                      
            if( activeStoreName.equals( "solr" ) ){
                // implicitly cleared by clearing index
            }

            if( activeStoreName.equals( "derby" ) ){
                recordStore = new DerbyRecordStore( psqContext );
                recordStore.clear();
            }
        } catch( MalformedURLException mux ){
            log.info( mux.getMessage(), mux );
        }              
    }

    //--------------------------------------------------------------------------

    public void processFile( String file ){
        File srcFl = new File( file );
        try{
            String cpath = srcFl.getCanonicalPath();
            this.root = cpath;
            enqueue( cpath, srcFl );        
        } catch( IOException ex ){ 
            log.info( ex.getMessage(), ex );
        }
    }
    
    public void processDirectory( String dir, boolean desc ){
        File dirFl = new File( dir );
        try{
            processFiles( dirFl.getCanonicalPath(), dirFl, desc );
        } catch( Exception ex ){
            log.info( ex.getMessage(), ex );
        }
    }
    
    private void processFiles( String root, File file, boolean desc ){
                
        this.root = root;
        File[] filesAndDirs = file.listFiles();
        List<File> filesDirs = Arrays.asList( filesAndDirs );
        
        if( ! file.isFile() ){
            log.info( "IndexBuilder: processing directory=" + file );
        }
        
        for( File sfile : filesDirs) {
            if ( ! sfile.isFile() && desc) {
                //recursive call: only when directory       
                processFiles( root, sfile, desc );
            } else {
                enqueue( root, sfile );                
            }
        }
    }
    
    List<List<File>> thq;

    int cth = 0;

    public void enqueue( String root, File file ){
        
        if ( file.isFile() ) {

            if(thq == null ){
                thq = new ArrayList<List<File>>();
            }       
            if( thq.size() <= cth ){
                thq.add( new ArrayList<File>());
            }
            thq.get( cth ).add( file );
            cth++;
            if( cth == builderTCount ) cth = 0;
        }
    }
    
    public void start(){
        
        // start threads
        
        List<IndexThread> itl = new ArrayList<IndexThread>();

        for( Iterator<List<File>> thi = thq.iterator(); thi.hasNext(); ){

            List<File> fq = thi.next();
            log.info( "IndexBuilder: queue size=" + fq.size() );
            for( Iterator<File> fqi = fq.iterator(); fqi.hasNext(); ){
                log.info( "IndexBuilder:  file: " + fqi.next() );
            }
            
            IndexThread it = new IndexThread( psqContext, host, solrconTCount,
                                              root, fq, format, zip );
            itl.add( it );

            log.info( "IndexBuilder:  starting thread..." );
            it.start();            
        }

        // wait till all threads are done
        //-------------------------------

        boolean running = true;

        try{
            while( running ){        

		running = false;

                for( Iterator<IndexThread> iti = itl.iterator(); 
                     iti.hasNext(); ){
                    IndexThread cit = iti.next();
                    if( cit.getState() != Thread.State.TERMINATED ){
                        running = true;
                    }
                    System.out.println( cit.getState() + ":" );
                }
                if( running ){
		    log.info( "IndexBuilder: still running...");
                    Thread.sleep( 1000*5 );
                } else {
		    log.info( "IndexBuilder: done...");
		}
            }

	    log.info( "IndexBuilder: all threads terminated");
            
        } catch( InterruptedException ix ){
            log.info( ix.getMessage(), ix );
        }
    }

    public void index( String root, File file ) throws IOException{
                
        String name = null;
        
        name = file.getCanonicalPath().replace( root, "" );
        
        log.info( "index: processing file:" + file );
        log.info( "index: file("
                  + "canonical=" + file.getCanonicalPath()
                  + "): " + name );
        log.info( "index: root:" + root );
        
        String compress = zip ? "zip" : "";
        
        // transform/add -> index
        //-----------------------

        if( recordIndex!= null ){
            recordIndex.addFile( file, name, format, compress );
        }
        // transform/add -> datastore
        //---------------------------
        if( recordStore!= null ){
            recordStore.addFile( file, name, format, compress );       
        }
    }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class IndexThread extends Thread{

    RecordIndex recordIndex = null;
    RecordStore recordStore = null;

    boolean zip = false;
    String root;
    String format;
    
    List<File> fileq;

    public IndexThread( PsqContext psqContext, String host, int solrconTCount,
                        String root, List<File> files, 
                        String format, boolean zip ){
        fileq = files;
        this.zip = zip;
        this.root = root;
        this.format = format;
        
        try{
            
            // SOLR Index
            //-----------
        
            String activeIndexName = psqContext.getActiveIndexName();
            
            if( activeIndexName.equals( "solr" ) ){
                recordIndex = new SolrRecordIndex( psqContext, host, 
                                                   solrconTCount );
                recordIndex.initialize();
                recordIndex.connect();
            }
            
            // Derby data store
            //-----------------
            
            String activeStoreName =  psqContext.getActiveStoreName();
            
            if( activeStoreName.equals( "derby" ) ){
                recordStore = new DerbyRecordStore( psqContext, host );
                recordStore.initialize();
            }
        } catch( MalformedURLException mux ){
            mux.printStackTrace();
        }
    }

    //--------------------------------------------------------------------------

    public void run(){
 
        Log log = LogFactory.getLog( this.getClass() );
        log.info( "IndexThread: processing fileq:" + fileq );

        String compress = zip ? "zip" : "";
        
        for( Iterator<File> fi = fileq.iterator(); fi.hasNext(); ){

            try{
                
                File file = fi.next();
                String name = file.getCanonicalPath().replace( root, "" );

                log.info( "IndexThread: processing file:" + file );
                log.info( "IndexThread: file("
                          + "canonical=" + file.getCanonicalPath() 
                          + "): " + name );
                log.info( "IndexThread: root:" + root );
                
                // transform/add -> index
                //-----------------------
                
                if( recordIndex != null){
                    recordIndex.addFile( file, name, format, compress );      
                }

                // transform/add -> datastore
                //---------------------------
                if( recordStore != null){  
                    recordStore.addFile( file, name, format, compress );
                }
            }catch( Exception ex ){
                log.info( ex.getMessage(), ex );
            }
        }
        
    }
}
