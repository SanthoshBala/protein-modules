package org.hupo.psi.mi.psicquic.server.store.derby;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # DerbyRecordStore: apache derby-backed RecordStore implementation
 #
 #
 #=========================================================================== */

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Clob;
import java.sql.PreparedStatement;

import org.w3c.dom.*;
import javax.xml.parsers.*;

import org.hupo.psi.mi.psicquic.server.*;
import org.hupo.psi.mi.psicquic.server.store.*;

import java.io.*;
import java.net.*;
import java.util.*;
import java.util.zip.*;

import java.util.regex.*;

import org.hupo.psi.mi.psicquic.util.JsonContext;
import org.hupo.psi.mi.psicquic.server.*;
import org.hupo.psi.mi.psicquic.server.store.*;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

//------------------------------------------------------------------------------

public class DerbyRecordStore implements RecordStore{

    public static final int TIMEOUT_LONG = 30;
    public static final int TIMEOUT_SHORT = 5;


    Log log;
    Connection dbcon = null;

    Map<String,Map<String,PsqTransformer>> inTransformerMap = null;

    String rmgrURL = null;

    PsqContext psqContext = null;
    String host = null;
    
    public DerbyRecordStore(){

        log = LogFactory.getLog( this.getClass() );

        try{
            Class.forName("org.apache.derby.jdbc.EmbeddedDriver");
        } catch( Exception ex ){
            ex.printStackTrace();
        }
    }

    public DerbyRecordStore( PsqContext context ){

        this.psqContext = context;
        log = LogFactory.getLog( this.getClass() );
        try{
            Class.forName("org.apache.derby.jdbc.EmbeddedDriver");
        } catch( Exception ex ){
            ex.printStackTrace();
        }
    }

    public DerbyRecordStore( PsqContext context, String host ){
        
        this.psqContext = context;
        log = LogFactory.getLog( this.getClass() );
        if( host != null ){
            this.host = host;
        }
        
        try{
            Class.forName("org.apache.derby.jdbc.EmbeddedDriver");
        } catch( Exception ex ){
            ex.printStackTrace();
        }
    }
    
    public void setPsqContext( PsqContext context ){
        this.psqContext = context;
    }
    
    public void initialize(){
        
    }
    
    private void connect(){
        
        if( dbcon == null ){
            
            Log log = LogFactory.getLog( this.getClass() );
            log.info( "DerbyRecordDao:connect" );
            
            if( psqContext != null && psqContext.getJsonConfig() != null ){

                Map derbyCfg = 
                    (Map) ((Map) psqContext.getJsonConfig().get("store"))
                    .get("derby");
                try{
                    String derbydb = (String) derbyCfg.get("derby-db");
                    log.info( "               location: " + derbydb );

                    dbcon = 
                        DriverManager.getConnection( "jdbc:derby:" + 
                                                     derbydb + ";create=true");
                } catch( Exception ex ){
                    ex.printStackTrace();
                }
                
                try{
                    Statement st = dbcon.createStatement();
                    st.setQueryTimeout( TIMEOUT_SHORT );
                    ResultSet rs = 
                        st.executeQuery( " select count(*) from record" );
                    while( rs.next() ){
                        log.info( "   record count= " + rs.getInt(1) );
                    }
                } catch( Exception ex ){
                    
                    // missing table ?
                    log.info( "   creating record table" ); 
                    create();
                }
            }
        }
    }

    //--------------------------------------------------------------------------

    private void create(){
        try{
            Statement st = dbcon.createStatement();
            st.setQueryTimeout( TIMEOUT_LONG ); 
            
            st.executeUpdate( "create table record " +
                              "(pk int generated always as identity," +
                              " rid varchar(256),"+
                              " format varchar(32), record clob )");
            
            st.executeUpdate( "create index r_rid on record (rid)" );
            
            st.executeUpdate( "create index r_ft on record (format)" );
            
        } catch( Exception ex ){
            ex.printStackTrace();
        }
    }
    
    public void shutdown(){
        if( dbcon != null ){      
            try{
                DriverManager.getConnection( "jdbc:derby:;shutdown=true" );
            } catch( Exception ex ){
                Log log = LogFactory.getLog( this.getClass() );
                log.info( "DerbyRecordDao(shutdown): " + ex);
            }
        }
    }
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    public void addRecord( String rid, String record, String format ){

        Log log = LogFactory.getLog( this.getClass() );
        log.debug("rid=" + rid + " format=" + format 
                  + " record(length)=" + record.length() );
        
        connect();
        try{
            PreparedStatement pst = dbcon
                .prepareStatement( "insert into record (rid, record, format)" +
                                   " values (?,?,?)" );
            
            pst.setString( 1, rid );
            
            if( record != null && format != null ){
                pst.setString( 2, record );
                pst.setString( 3, format );
            
                pst.executeUpdate();
            } 
        }catch( Exception ex ){
            ex.printStackTrace();
        }
    }

    //--------------------------------------------------------------------------

    public String getRecord( String rid, String format ){
        connect();
        String record = "";
        
	Log log = LogFactory.getLog( this.getClass() );

        try{
            PreparedStatement pst = dbcon
                .prepareStatement( "select rid, record, format from record" +
                                   " where rid = ? and format= ?" );
            
            pst.setString( 1, rid );
            pst.setString( 2, format );
            ResultSet rs =  pst.executeQuery();

            String rt = "";

            while( rs.next() ){
                Clob rc = rs.getClob( 2 );
                record = rc.getSubString( 1L, 
                                          new Long(rc.length()).intValue() );
                rt = rs.getString(3);
            } 

	log.info( "DerbyRecordDao(getRecord): recId=" + rid + " rt=" + rt 
                  + "  record=" + record );
        
        }catch( Exception ex ){
            ex.printStackTrace();
        }
        return record;
    }

    //--------------------------------------------------------------------------
    
    public List<String> getRecordList( List<String> rid, String format ){
        
        connect();
        List<String> recordList = new ArrayList<String>();
        
        try{
            PreparedStatement pst = dbcon
                .prepareStatement( "select rid, record from record" +
                                   " where rid = ? and format = ?" );
            
            for( Iterator<String> i = rid.iterator(); i.hasNext(); ){
                pst.setString( 1, i.next() );
                pst.setString( 2, format );
                ResultSet rs =  pst.executeQuery();
                
                while( rs.next() ){
                    Clob rc = rs.getClob(2);
                    String record = 
                        rc.getSubString( 1L, 
                                         new Long(rc.length()).intValue() );
                    recordList.add( record );
                }
            }
        } catch( Exception ex ){
            ex.printStackTrace();
        }
        return recordList;
    }


    //--------------------------------------------------------------------------

    public String toString( org.hupo.psi.mi.psicquic.server.ResultSet rset ){
        
        String rstr = "";
        
        return rstr;
    }
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    public void addFile( File file, String fileName, 
                         String format, String compress){
                
        Map trCfg = (Map) ((Map) psqContext.getJsonConfig().get("store"))
            .get("transform");
        
        List rtrList = (List) trCfg.get( format );

        log.debug( "DerbyRecordStore:" 
                   + " addFile: file=" + file + " name=" + fileName );
        log.debug( "DerbyRecordStore:"
                   + " addFile: format=" + format + " trl=" + rtrList );
        
        if( inTransformerMap == null ){
            inTransformerMap 
                = new HashMap<String,Map<String,PsqTransformer>>();            
        }
        
        Map<String,PsqTransformer> itm = inTransformerMap.get( format );
        if( itm == null ){
            itm = new HashMap<String,PsqTransformer>();
            inTransformerMap.put( format, itm );
            itm = inTransformerMap.get( format );
        }
        
        if( rtrList != null ){
            for( Iterator it = rtrList.iterator(); it.hasNext(); ){
                
                Map itr = (Map) it.next();
                log.debug( "DerbyRecordStore:" +
                           " addFile: view=" + itr.get( "view" ) );
                if( ((String) itr.get("type")).equalsIgnoreCase("XSLT") &&
                    (Boolean) itr.get("active") ){

                    if( itm.get( itr.get("view") ) == null ){            
                        
                    // initialize transformer
                    //-----------------------
                        
                        log.info( " Initializing transformer: format=" + format
                                  + " type=XSLT config=" + itr.get("config") );
                        
                        PsqTransformer rt =
                            new XsltTransformer( (Map) itr.get("config") );
                        
                        itm.put( (String) itr.get("view"), rt );
                    }

                    PsqTransformer rt = itm.get( itr.get("view") );                     
                    
                    try{
                        if( compress!= null 
                            && compress.equalsIgnoreCase("zip") ){
                            
                            processZipFile( rt, (String) itr.get( "view" ), 
                                            fileName, new ZipFile( file ));                            
                        } else {
                            processFile( rt, (String) itr.get( "view" ),
                                         fileName, new FileInputStream( file ));
                        }
                        
                    } catch( Exception ex ){
                        log.info( ex.getMessage(), ex );
                        return;
                    }
                }                
            }            
        }
    }

    //--------------------------------------------------------------------------

    private void processZipFile( PsqTransformer rt, String viewName,
                                 String fileName,  ZipFile zf )
        throws java.io.IOException{
        
        for( Enumeration zfe = zf.entries(); zfe.hasMoreElements(); ){
            
            ZipEntry ze = (ZipEntry) zfe.nextElement();
            if( !ze.isDirectory() ){   
                InputStream is = zf.getInputStream( ze );
                processFile( rt, viewName,
                             fileName + "::" + ze.getName() , is );
            }
        }
    }

    //--------------------------------------------------------------------------

    private void processFile( PsqTransformer rt, String viewName, 
                              String fileName, InputStream is ){
        
        rt.start( fileName, is );

        while( rt.hasNext() ){
					
            Map cdoc= rt.next();
            String recId = (String) cdoc.get( "recId" );
            NodeList fl = (NodeList) cdoc.get( "dom" );
                            
            String vStr = null;
            
            for( int j = 0; j< fl.getLength() ;j++ ){
                if( fl.item(j).getNodeName().equals( "field") ){
                    Element fe = (Element) fl.item(j);
                    String name = fe.getAttribute("name");
                    String value = fe.getFirstChild().getNodeValue();
                    
                    log.debug( " Node: name=" + name 
                               + " length(value)="+ value.length() );
                    
                    if( name.equals( "recId" ) ){
                        recId = value;
                    }
                    if( name.equals( "view" ) ){
                        vStr = value;
                    }
                }
            }
            
            try{
                this.add( recId, viewName, vStr );
            }catch( Exception ex ){
                ex.printStackTrace();
            }
        }
    }

    //--------------------------------------------------------------------------
    
    private void add( String pid, String vType, String view ){
        
        log.info( "PID=" + pid ); 
        log.info( "  VTP=" + vType ); 
        log.info( "  VIEW=" + view.substring(0,24) + "..." ); 

        // add record
        //-----------
        
        try{
            String postData = URLEncoder.encode("op", "UTF-8") + "="
                + URLEncoder.encode( "add", "UTF-8");
            postData += "&" + URLEncoder.encode("pid", "UTF-8") + "="
                + URLEncoder.encode( pid, "UTF-8");
            postData += "&" + URLEncoder.encode("vt", "UTF-8") + "="
                + URLEncoder.encode( vType, "UTF-8");
            
            postData += "&" + URLEncoder.encode("vv", "UTF-8") + "="
                + URLEncoder.encode( view, "UTF-8");

            if( rmgrURL == null ){
                rmgrURL = (String) 
                    ((Map) psqContext.getJsonConfig().get("store"))
                    .get("record-mgr");            

                if( host != null ){
                    rmgrURL = hostReset( rmgrURL, host );
                }
            }
            
            URL url = new URL( rmgrURL );
            URLConnection conn = url.openConnection();
            conn.setDoOutput( true );
            OutputStreamWriter wr =
                new OutputStreamWriter( conn.getOutputStream() );
            wr.write(postData);
            wr.flush();
            
            // Get the response

            InputStreamReader isr = 
                new InputStreamReader( conn.getInputStream() );
            BufferedReader rd = new BufferedReader( isr );
            String line = null;

            log.info( "  Response:" );
            while ((line = rd.readLine()) != null) {           
                log.info( line );
            }
            wr.close();
            rd.close();
        } catch (Exception e) {
        }
    }

    //--------------------------------------------------------------------------

    public void clear(){

        log = LogFactory.getLog( this.getClass() );
        
        try{
            String postData = URLEncoder.encode("op", "UTF-8") + "="
                + URLEncoder.encode( "clear", "UTF-8");
            
            if( rmgrURL == null ){
                rmgrURL = (String) 
                    ((Map) psqContext.getJsonConfig().get("store"))
                    .get("record-mgr");            

                if( host != null ){
                    rmgrURL = hostReset( rmgrURL, host );
                }
            }
            
            URL url = new URL( rmgrURL );
            URLConnection conn = url.openConnection();
            conn.setDoOutput( true );
            OutputStreamWriter wr =
                new OutputStreamWriter( conn.getOutputStream() );
            wr.write(postData);
            wr.flush();
            
            // Get the response
            BufferedReader rd =
                new BufferedReader( new InputStreamReader( conn
                                                           .getInputStream()));
            String line;
            
            log.info( "  Response:" );
            while ((line = rd.readLine()) != null) {           
                log.info( line );
            }
            wr.close();
            rd.close();
        } catch (Exception e) {
        }
    }

    //--------------------------------------------------------------------------
    
    public void clearLocal(){
        
        log = LogFactory.getLog( this.getClass() );
        
        if( dbcon == null ){
            connect();
        }

        try{
            Statement st = dbcon.createStatement();
            st.setQueryTimeout(60);
            st.executeUpdate( "truncate table record");
            log.info( "record table truncated" );
        } catch( Exception ex ){
            // missing table ?
            log.info( ex.getMessage(), ex );
            log.info( "creating record table" );
            create();
        }
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
                return newUrl;

            } else {
                return url;
            }
        } catch(Exception ex ){
            return url;
        }
    }
}
