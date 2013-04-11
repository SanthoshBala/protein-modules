package org.hupo.psi.mi.psicquic.util;

/* =============================================================================
 # $Id:: JsonContext.java 937 2012-05-29 15:39:09Z lukasz99                    $
 # Version: $Rev:: 937                                                         $
 #==============================================================================
 #                                                                             $
 # JsonContext: JSON-based configuration                                       $
 #                                                                             $
 #     TO DO:                                                                  $
 #                                                                             $
 #=========================================================================== */

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.*;
import java.io.*;
import org.json.*;

public class JsonContext {

    public JsonContext() {}

    Map<String,Object> config;

    private static final int BUFFER_SIZE = 4096;

    //--------------------------------------------------------------------------

    public void setConfig( Map<String,Object> config ) {
	this.config = config;
    }

    public Map<String,Object> getConfig() {
	return config;
    }

    //--------------------------------------------------------------------------
    // JSON-based configuration
    //-------------------------

    private String jsonConfigString;

    public String getJsonConfigString() {
        return jsonConfigString;
    }

    public void setJsonConfigString( String def ) {
        this.jsonConfigString = def;
    }

    private JSONObject jsonConfigObject;

    public JSONObject getJsonConfigObject() {
        return jsonConfigObject;
    }

    public void setJsonConfigObject( JSONObject def ) {
        this.jsonConfigObject = def;
    }

    private Map<String,Object> jsonConfigUtil = null;

    public Map<String,Object>  getJsonConfig() {

        if( jsonConfigUtil == null ){
            try{
                readJsonConfigDef();
            }catch( IOException ex ){
                Log log = LogFactory.getLog( JsonContext.class );
                log.info("JsonContext: configuration error");
            }
        }
        return jsonConfigUtil;
    }

    //--------------------------------------------------------------------------
    // read JSON definition
    //---------------------
   
    public void readJsonConfigDef() throws IOException{

        String jsonPath = null;

        if( config != null && config.get( "json-config" ) != null ){
            jsonPath = (String) config.get( "json-config" );
        }
        
	File jsfFile = new File( jsonPath );
	InputStream jsfIStr = null;

	if( !jsfFile.canRead() ){
	    jsfIStr = this.getClass().getClassLoader()
		.getResourceAsStream( jsonPath );
	} else {
	    jsfIStr = new FileInputStream( jsonPath );
	}

        //InputStream is = 
        //    ApplicationContextProvider.getResourceAsStream( jsonPath );
        readJsonConfigDef( jsfIStr );
    }
   

    public void readJsonConfigDef( InputStream is ) {

	Log log = LogFactory.getLog( JsonContext.class );

	StringBuffer sb = new StringBuffer();
	char[] buffer = new char[BUFFER_SIZE];

	try {
	    InputStreamReader ir = new InputStreamReader( is );
	    int len =0;
	    while ( (len = ir.read( buffer, 0, BUFFER_SIZE ) ) >= 0 ) {
		sb.append( buffer , 0, len);
	    }

	} catch ( Exception e ) {
	    log.info( e.toString() );
	}

	String jsonConfigDef = sb.toString();
	log.info( "unparsed=" + jsonConfigDef);
	try {
	    JSONObject jo = new JSONObject( jsonConfigDef );
	    log.info( "parsed: " +jo.toString() );
	    jsonConfigObject = jo;
	    jsonConfigString = jo.toString();
	    jsonConfigUtil = json2util( jo );
	} catch ( JSONException jex ) {
	    log.info( "parsing error: " + jex.toString() );
	}
    }


    //--------------------------------------------------------------------------
    // write JSON definition
    //----------------------

    public void writeJsonConfigDef( PrintWriter pw ) {
		Log log = LogFactory.getLog( this.getClass() );
		JSONObject currentJSONConfigObject = util2json(jsonConfigUtil);
		try {
			pw.print(currentJSONConfigObject.toString(2));
		} catch (JSONException jex) {
	    	log.info( "JSON printing error: " + jex.toString());
		}
    }

    //--------------------------------------------------------------------------
    // java.util -> JSON conversion
    //-----------------------------

    public JSONObject util2json( Map<String,Object> util ) {

	if ( util != null ) {
	    return new JSONObject( util );
	} else {
	    return new JSONObject();
	}
    }

    //--------------------------------------------------------------------------
    // JSON -> java.util conversion
    //-----------------------------

    public Map<String,Object> json2util( JSONObject json ) {

	Map<String,Object> util = new HashMap<String,Object>();

	for ( Iterator ki = json.keys(); ki.hasNext(); ) {
	    String key = (String) ki.next();
	    Object o = null;
	    try {
		o = json.get( key );
	    } catch ( JSONException jex ) {
		key = null;
	    }
	    if ( key != null &&
                 o.getClass().getName().equals( "org.json.JSONObject" ) ) {

                // JSONObject
                o = json2util( (JSONObject) o );
	    }

	    if ( key != null &&
                 o.getClass().getName().equals( "org.json.JSONArray" ) ) {

                // JSONArray
		o = json2util( (JSONArray) o );
	    }
	    if (key != null ){
		util.put( key,o );
	    }
	}
	return util;
    }

    //--------------------------------------------------------------------------

    public List json2util( JSONArray json ) {

	List util = new ArrayList();

	for ( int i = 0; i < json.length(); i++ ) {

	    Object o = null;
	    try {
	        o = json.get( i );
	    } catch ( JSONException jex ) {
		o = null;
	    }

	    if ( o != null &&
                 o.getClass().getName().equals( "org.json.JSONObject" ) ) {

                // JSONObject
                o = json2util( (JSONObject) o );
	    }

	    if ( o != null &&
                 o.getClass().getName().equals( "org.json.JSONArray" ) ) {

                // JSONArray
		o = json2util( (JSONArray) o );
	    }
	    if (o != null ){
		util.add( o );
	    }
	}
	return util;
    }
}
