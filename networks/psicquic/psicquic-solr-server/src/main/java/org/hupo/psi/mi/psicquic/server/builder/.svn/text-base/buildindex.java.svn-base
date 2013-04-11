package org.hupo.psi.mi.psicquic.server.builder;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # BuildIndex: Build index
 #
 #=========================================================================== */

import java.util.*;
import java.net.*;
import java.io.*;

import java.lang.Thread.State;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.apache.commons.cli.*;

//import org.w3c.dom.*;
//import javax.xml.parsers.*;

import org.hupo.psi.mi.psicquic.util.JsonContext;

public class buildindex{
    
    public static final String 
        DEFAULT_CONTEXT = "/etc/psq-context-default.json";
    
    public static final String 
        DEFAULT_RECORD_FORMAT = "mif254";   

    public static String rfrmt = DEFAULT_RECORD_FORMAT;
    
    public static String host = null;
    public static String idir = null;
    public static String ifile = null;
    
    public static boolean clr = false;
    public static boolean zip = false;
    public static boolean desc = false;

    public static int btCount = -1;
    public static int stCount = -1;
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    public static void main( String [ ] args ){
        
        Options options = new Options();
        
        Option clrOption = OptionBuilder.withLongOpt( "clear" )
            .withDescription( "clear index" )
            .create( "clr" );

        options.addOption( clrOption );

        Option ctxOption = OptionBuilder.withLongOpt( "context" )
            .withArgName( "file.json" ).hasArg()
            .withDescription( "configuration file" )
            .create( "ctx" );

        options.addOption( ctxOption );
        
        Option dirOption = OptionBuilder.withLongOpt( "dir" )
            .withArgName( "dir" ).hasArg()
            .withDescription( "input file directory" )
            .create( "d" );
        
        options.addOption( dirOption );

        Option fileOption = OptionBuilder.withLongOpt( "file" )
            .withArgName( "file" ).hasArg()
            .withDescription( "input file" )
            .create( "f" );

        options.addOption( fileOption );
        
        Option fmtOption = OptionBuilder.withLongOpt( "format" )
            .withArgName( "format" ).hasArg()
            .withDescription( "input record format" )
            .create( "fmt" );
        
        options.addOption( fmtOption );
        
        Option hlpOption = OptionBuilder.withLongOpt( "help" )
            .withDescription( "help " )
            .create( "help" );
        
        options.addOption( hlpOption );

        Option hostOption = OptionBuilder.withLongOpt( "host" )
            .withArgName( "host" ).hasArg()
            .withDescription( "host" )
            .create( "h" );

        options.addOption( hostOption );
               
        Option recOption = OptionBuilder.withLongOpt( "r" )
            .withDescription( "recursively process directory" )
            .create( "r" );
        
        options.addOption( recOption );

        Option tbOption = OptionBuilder.withLongOpt( "threads-builder" )
            .withArgName( "count" ).hasArg()
            .withDescription( "number of builder threads" )
            .create( "tb" );
        
        options.addOption( tbOption );
        
        Option tsOption = OptionBuilder.withLongOpt( "threads-solr" )
            .withArgName( "count" ).hasArg()
            .withDescription( "number of solr connection threads" )
            .create( "ts" );
        
        options.addOption( tsOption );
        
        Option zipOption = OptionBuilder.withLongOpt( "zip" )
            .withDescription( "zipped files" )
            .create( "z" );

        options.addOption( zipOption );
        
        String context = null;     
        String ifrmt = DEFAULT_RECORD_FORMAT;
        
        try{
            CommandLineParser parser = new PosixParser();
            CommandLine cmd = parser.parse( options, args);

            if( cmd.hasOption("help") ){

                HelpFormatter formatter = new HelpFormatter();
                formatter.setWidth( 127 );
                formatter.printHelp( "BuildSolrDerbyIndex", options );
                System.exit(0);
            }
            
            if( cmd.hasOption( "clr" ) ){
                clr = true;
            } 
            
            if( cmd.hasOption( "ctx" ) ){
                context = cmd.getOptionValue( "ctx" );
            } 
        
            if( cmd.hasOption( "fmt" ) ){
                ifrmt = cmd.getOptionValue( "fmt" );
            }

            if( cmd.hasOption( "d" ) ){
                idir = cmd.getOptionValue( "d" );
            }
            
            if( cmd.hasOption( "f" ) ){
                ifile = cmd.getOptionValue( "f" );
            }
            
            if( cmd.hasOption( "h" ) ){
                host = cmd.getOptionValue( "h" );
            }
            
            if( cmd.hasOption( "z" ) ){
                zip = true;
            }
            
            if( cmd.hasOption("r") ){
                desc = true;
            }

            if( cmd.hasOption("tb") ){
                try{
                    String sBtc =  cmd.getOptionValue( "tb" );
                    int btc = Integer.parseInt( sBtc ); 
                    btCount = btc;
                }catch( NumberFormatException nfx ){
                    System.out.println( "TB format Error. Using default" );
                }
            }

            if( cmd.hasOption("ts") ){
                try{
                    String sStc =  cmd.getOptionValue( "ts" );
                    int stc = Integer.parseInt( sStc ); 
                    stCount = stc;
                }catch( NumberFormatException nfx ){
                    System.out.println( "TS format Error. Using default" );
                }
            }



        } catch( Exception exp ) {
            System.out.println( "BuildIndex: Options parsing failed. " +
                                "Reason: " + exp.getMessage() );
            HelpFormatter formatter = new HelpFormatter();
            formatter.setWidth(127);
            formatter.printHelp( "BuildIndex", options );
            System.exit(1);
        }

        IndexBuilder ibuilder = null;
        
        if( context != null ){
            System.out.println( "Context: " + context );
            ibuilder = new IndexBuilder( context, host, btCount, stCount, 
                                         ifrmt, zip );
        }else{
            
            System.out.println( "Context(default): " + DEFAULT_CONTEXT );
            try{
                InputStream ctxStream = buildindex.class
                    .getResourceAsStream( DEFAULT_CONTEXT );
                
                ibuilder = new IndexBuilder( ctxStream, host, btCount, stCount, 
                                             ifrmt, zip );
            } catch( Exception ex){
                ex.printStackTrace();
            }
            
        }

        if( clr ){
            System.out.println( "Clear index" );
            ibuilder.clear();
        }

        if( ibuilder != null && ifile != null ){  
            ibuilder.processFile( ifile );
            ibuilder.start();
        } else { // file option has a precedence over directory
            if( ibuilder != null && idir != null ){
                ibuilder.processDirectory( idir, desc );
                ibuilder.start();
            }
        }
    }
}
