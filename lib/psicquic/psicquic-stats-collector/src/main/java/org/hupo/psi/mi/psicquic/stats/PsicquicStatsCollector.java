package org.hupo.psi.mi.psicquic.stats;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.gdata.client.spreadsheet.SpreadsheetService;
import com.google.gdata.data.spreadsheet.Cell;
import com.google.gdata.data.spreadsheet.CellEntry;
import com.google.gdata.data.spreadsheet.SpreadsheetEntry;
import com.google.gdata.util.ServiceException;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.exception.ExceptionUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.registry.client.PsicquicRegistryClientException;
import org.hupo.psi.mi.psicquic.registry.client.registry.DefaultPsicquicRegistryClient;
import org.hupo.psi.mi.psicquic.registry.client.registry.PsicquicRegistryClient;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;
import org.springframework.beans.factory.BeanFactory;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import org.springframework.mail.MailSender;
import org.springframework.mail.SimpleMailMessage;
import psidev.psi.mi.tab.PsimiTabException;
import psidev.psi.mi.tab.PsimiTabReader;
import psidev.psi.mi.tab.model.BinaryInteraction;
import psidev.psi.mi.tab.model.CrossReference;
import uk.ac.ebi.intact.google.spreadsheet.SpreadsheetFacade;
import uk.ac.ebi.intact.google.spreadsheet.WorksheetFacade;

import java.io.*;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.*;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// How to install the Gdata libs locally ?
//
// 1. Download the gdata libs from Google Code
//
// 2. unpack and install the following jar files present in the gdata/java/libs directory:
//    gdata-spreadsheet-3.0.jar
//    gdata-core-1.0.jar
//
//    Command lines:
//
//        mvn install:install-file -Dfile=gdata-spreadsheet-3.0.jar -DgroupId=com.google.gdata.client.spreadsheet \
//                                 -DartifactId=gdata-spreadsheet -Dversion=3.0 -Dpackaging=jar
//
//        mvn install:install-file -Dfile=gdata-core-1.0.jar -DgroupId=com.google.gdata.data \
//                                 -DartifactId=gdata-core -Dversion=1.0 -Dpackaging=jar
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Utility that collect PSICQUIC statistics from existing services and publishes it into a Google Spreadsheet.
 * <p/>
 * The tool can be used to monitor the PSICQUIC services hosted on on the main Registry but an other Registry can also
 * be used by passing the System property 'psicquic.registry.url'. Whether it's a totally different registry or the
 * public one with some exclusion rules or tag filtering.
 * <p/>
 * One can use the default SMTP configuration file provided (resources/META-INF/smtp.properties) or use an external one
 * by setting it's location as a System's property: 'smtp.config.file'.
 * <p/>
 * Note: We used a lot of code sample from:
 * http://code.google.com/apis/spreadsheets/data/3.0/developers_guide_java.html
 * http://ga-api-java-samples.googlecode.com/svn/trunk/src/v1/SpreadsheetExporter.java
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class PsicquicStatsCollector {

    private static final Log log = LogFactory.getLog( PsicquicStatsCollector.class );

    public static final String FILE_SEPARATOR = System.getProperty( "file.separator" );

    private static SimpleDateFormat DATE_FORMAT = new SimpleDateFormat( "yyyy-MM-dd" );
    private static SimpleDateFormat DATE_AS_DIR = new SimpleDateFormat( "yyyy/MM/dd" );

    private static final String PSICQUIC_REGISTRY_URL_KEY = "psicquic.registry.url";

    private static final String SMTP_CONFIG_FILE_KEY = "smtp.config.file";

    private static final int PSICQUIC_DEFAULT_TIMEOUT = 60000;
    private static final int PSICQUIC_BATCH_SIZE = 500;

    private static final String TOTAL_COLUMN = "Total";
    private static final String DATE_COLUMN = "Date";

    private MailSender mailSender;
    private String senderEmail;
    private String mailSubjectPrefix;

    private List<String> recipients;
    private Config config;

    private int numberOfTests = 10;

    public PsicquicStatsCollector() throws IOException {
        // Initialize Spring for emails
        log.info("Initializing Spring ...");
        String[] configFiles = new String[]{"/META-INF/beans.spring.xml"};
        BeanFactory factory = new ClassPathXmlApplicationContext( configFiles );
        log.info("Spring was initialized.");

        // inject my spring beans
        this.mailSender = ( MailSender ) factory.getBean( "mailSender" );
        this.config = ( Config ) factory.getBean( "statsConfig" );

        log.info("Config post spring initialization:" + config);

        // load email properties
        File smtpConfig = null;
        if( ! config.hasSmtpConfigFile() ) {
            log.info( "Using default SMTP config from classpath: '"+ config.DEFAULT_SMTP_CONFIG +"'" );
            smtpConfig = new File( PsicquicStatsCollector.class.getResource( config.DEFAULT_SMTP_CONFIG ).getFile() );
        } else {
            smtpConfig = new File( config.getSmtpConfigFile() );
        }
        log.info( "Using SMTP config from: '"+ smtpConfig.getAbsolutePath() +"'" );

        Properties properties = new Properties();
        FileInputStream in = new FileInputStream( smtpConfig );
        properties.load( in );
        in.close();

        mailSubjectPrefix = properties.getProperty( "email.subject.prefix" );
        log.info( "mailSubjectPrefix = " + mailSubjectPrefix );
        String recipientList = properties.getProperty( "email.recipients" );

        recipients = new ArrayList( Arrays.asList( recipientList.split( "," ) ) );
        log.info( "recipients = " + recipients );

        senderEmail = properties.getProperty( "email.sender" );
        log.info( "senderEmail = " + senderEmail );
        if( StringUtils.isEmpty(senderEmail ) ) {
            // disable email sending facilities
            log.warn( "No email sender specified, running with email facilities disabled." );
            senderEmail = null;
            mailSender = null;
        }

        log.info( "configuration post initialization: " + config );
    }

    ///////////////////////////////////////
    // Mail handling

    public void sendEmail( String title, String body ) {
        if ( mailSender == null ) {
            log.debug( "------------------------------------------------------------" );
            log.debug( "From:  " + senderEmail );
            log.debug( "To:    " + recipients );
            log.debug( "Title: " + mailSubjectPrefix + title );
            log.debug( body );
            log.debug( "------------------------------------------------------------" );
        } else {
            final SimpleMailMessage message = new SimpleMailMessage();
            message.setTo( recipients.toArray( new String[]{} ) );
            message.setFrom( senderEmail );
            message.setSubject( mailSubjectPrefix + " " + title );
            message.setText( body );
            mailSender.send( message );
        }
    }

    ///////////////////////////////////////
    // Google Spreadsheet related methods

    public List<String> updateSheet( final SpreadsheetService service,
                                     SpreadsheetEntry spreadsheetEntry,
                                     String worksheetName,
                                     Map<String, Long> db2count ) throws IOException, ServiceException {

        final List<String> updatedServices = Lists.newArrayList();

        final WorksheetFacade myWorksheet = WorksheetFacade.getWorksheetByName( service, spreadsheetEntry, worksheetName );
        if( myWorksheet == null ) {
            final String msg = "Could not find the worksheet '" + worksheetName + "' in spreadsheet '" + spreadsheetEntry.getTitle().getPlainText() + "'";
            sendEmail(msg, "");
            throw new IllegalArgumentException( msg );
        }

        Map<String, Integer> columnName2Index = myWorksheet.getColumnNameIndex();
        final int nextEmptyRow = myWorksheet.getNextEmptyRow();
        log.info( "Next empty Row: " + nextEmptyRow );
        final List<CellEntry> updatedCells = Lists.newArrayList();

        // If a total column is present, then sum up all cells into it
        boolean hasTotalColumn = false;
        int totalColumnIndex = -1;
        int currentTotalSum = 0;

        for ( Map.Entry<String, Integer> e : columnName2Index.entrySet() ) {
            final String db = e.getKey();
            final Integer colIndex = e.getValue();

            if ( colIndex == null ) {
                throw new IllegalArgumentException( "You must give a non null colIndex" );
            }

            final Long dbCount = db2count.get( db );

            log.info( "Updating worksheet for '"+ db +"' having count of "+ dbCount );

            String count;
            if( dbCount == null ) {

                if( TOTAL_COLUMN.equalsIgnoreCase( db ) ) {
                    totalColumnIndex = colIndex;
                    hasTotalColumn = true;

                    log.info( "Skipping column: " + db );
                    continue; // go to next databse
                }

                if( DATE_COLUMN.equalsIgnoreCase( db ) ) {
                    log.info( "Skipping column: " + db );
                    continue; // go to next databse
                }
            }

            CellEntry previousCell = null;
            CellEntry previousDateCell = null;

            // TODO create WorksheetFacade.hasData() that returns true if there is at least one row of data in it.
            if ( nextEmptyRow > 2 ) {
                previousCell = myWorksheet.getCell( nextEmptyRow - 1, colIndex );
                log.info( "Collecting data for previous cell (value="+ previousCell.getCell().getValue() +")" );
                previousDateCell = myWorksheet.getCell( nextEmptyRow - 1, 1 );
            }

            CellEntry cell = null;

            if ( dbCount != null ) {

                count = dbCount.toString();

                // we have data for that database
                cell = myWorksheet.getCell( nextEmptyRow, colIndex );

                if ( previousCell != null ) {

                    final String previousValue = previousCell.getCell().getValue();

                    if( previousValue == null ) {
                        log.info( "No previous cell available." );
                        if ( ! "0".equals( count ) ) {
                            if ( log.isDebugEnabled() ) log.info( worksheetName + ": statistics for " +
                                    db + " has changed since " +
                                    previousDateCell.getCell().getValue() );
                            // the service has a positive count so we record it in the stats
                            cell.changeInputValueLocal( count );
                            updatedServices.add( db );
                        } else {
                            log.info( "Skip update as the count is 0" );
                        }
                    } else {
                        log.info( "Previous cell available: " + previousValue );
                        if( ! previousValue.equals( count ) && ! "0".equals( count ) ) {
                            // the service has a positive count so we record it in the stats
                            cell.changeInputValueLocal( count );
                            updatedServices.add( db );
                        } else {
                            log.info( "Skip update as the count is 0 or the same as previous ("+ previousValue +")" );
                            // we store the value anyways but do not register the service as updated.
                            cell.changeInputValueLocal( count );
                        }
                    }

                } else {
                    // we don't have data yet so all services have to be updated
                    log.info( "No previous cell available." );
                    cell.changeInputValueLocal( count );
                    updatedServices.add( db );
                }

                if( ! StringUtils.isEmpty( cell.getCell().getInputValue() ) ) {
                    log.info( "Adding new cell for service '"+ db +"' with value: " + cell.getCell().getInputValue() );
                    updatedCells.add( cell );
                }

            } else {

                log.info( "++++ No data for " + db );
                log.info( "++++ previousCell: " + previousCell.getCell().getValue() );

                // We haven't got data for that database, copy the previous cell's data

                if ( previousCell != null ) {

                    final String previousValue = previousCell.getCell().getValue();
                    if( ! "0".equals( previousValue ) ) {
                        if ( log.isDebugEnabled() ) log.info( worksheetName + ": statistics for " + db +
                                " has NOT changed since " +
                                previousDateCell.getCell().getValue() +
                                " copying previous one instead: " + previousValue );

                        cell = myWorksheet.getCell( nextEmptyRow, colIndex.intValue() );
                        cell.changeInputValueLocal( previousCell.getCell().getValue() );
                        updatedCells.add( cell );
                    }

                } else {
                    log.info( "Previous cell is null :(" );
                }
            }

            // we may ommit part of the totalSum as we only count those that have new data.
            if( cell != null && ! StringUtils.isEmpty( cell.getCell().getInputValue() ) ) {
                log.info( "Adding " + cell.getCell().getInputValue() + " to total sum" );
                // keep track of current total
                currentTotalSum += Integer.parseInt( cell.getCell().getInputValue() );
            } else {
                if( previousCell != null ) {
                    log.info( "Adding " + previousCell.getCell().getInputValue() + " to total sum (previous cell)" );
                    // keep track of current total
                    if(previousCell.getCell().getValue() != null ) {
                        currentTotalSum += Integer.parseInt( previousCell.getCell().getValue() );
                    }
                }
            }

        } // All databases registered in the worksheet header

        if( updatedServices.isEmpty() ) {
            // no services updated, flush the collection of updatedCells
            updatedCells.clear();
        }

        final boolean hasNewData = ! updatedCells.isEmpty();

        if( hasTotalColumn ) {
            // process all existing rows and update the total column

            log.info( "Updating Total..." );
            log.info( "getNextEmptyRow(): " + myWorksheet.getNextEmptyRow() );
            for ( int row = 2; row < myWorksheet.getNextEmptyRow(); row++ ) {

                int totalSum = 0;
                boolean abortRow = false;

                log.info( "getNextEmptyColumn(): " + myWorksheet.getNextEmptyColumn() );
                for ( int col = 2; col < myWorksheet.getNextEmptyColumn(); col++ ) {
                    if( col != totalColumnIndex ) {
                        final CellEntry cell = myWorksheet.getCell( row, col );
                        final String value = cell.getCell().getValue();
                        if( value != null ) {
                            try {
                                final int valueInt = Integer.parseInt( value );
                                totalSum += valueInt;

                                log.info( myWorksheet.getCell( row, 1 ).getCell().getValue() + " -> " + valueInt);
                            } catch ( NumberFormatException e ) {
                                // report cell format issue by email, and skip
                                sendEmail( "Cell format error in worksheet: " + worksheetName,
                                        "Cell["+row+","+col+"]='"+value + "' when an integer was expected" );
                                abortRow = true;
                                continue; // exit column loop
                            }
                        }
                    } else {
                        log.info( "Total contains: " + myWorksheet.getCell( row, col ).getCell().getValue() );
                    }
                } // columns

                log.info( myWorksheet.getCell( row, 1 ).getCell().getValue() + " -> " + totalSum + " (SUM)");

                if( ! abortRow ) {
                    final CellEntry cell = myWorksheet.getCell( row, totalColumnIndex );
                    cell.changeInputValueLocal( String.valueOf( totalSum ) );
                    updatedCells.add( cell );
                }
            } // rows

            log.info( "currentTotalSum = " + currentTotalSum );
            if( hasNewData && currentTotalSum != 0 ) {
                // given that the last row of this sheet is not available in the WorksheetFacade until we save the updatedCells,
                // we save the last total additionally
                log.info( "Current total: " + currentTotalSum );
                final CellEntry cell = myWorksheet.getCell( myWorksheet.getNextEmptyRow(), totalColumnIndex );
                cell.changeInputValueLocal( String.valueOf( currentTotalSum ) );
                updatedCells.add( cell );
            }
        } // hasTotalColumn


        // Check if a database is not yet registered in the spreadsheet.
        for ( Map.Entry<String, Long> e : db2count.entrySet() ) {
            String db = e.getKey();

            final Integer colIndex = columnName2Index.get( db );
            if( colIndex == null ) {
                sendEmail( "Missing DB name in header of worksheet '" + worksheetName + "'",
                        "Please add " + db + " in the header of the spreadsheet - please visit "
                                + spreadsheetEntry.getSpreadsheetLink().getHref() );
            }
        }

        if ( log.isDebugEnabled() ) {
            log.info( "has new data ? " + hasNewData );
            log.info( "About to update spreadsheet with " + updatedCells.size() + " cell(s)" );
            for ( CellEntry c : updatedCells ) {
                final Cell cc = c.getCell();
                log.info( "["+cc.getCol() +","+cc.getRow()+"]: "+cc.getInputValue()+( cc.getValue() != null ? " (was "+ cc.getValue() +")" : " (new)" ) );
            }
        }

        if ( ! updatedCells.isEmpty() ) {
            if( hasNewData ) {
                // add date to the row
                final CellEntry cell = myWorksheet.getCell( nextEmptyRow, 1 );
                cell.changeInputValueLocal( today() );
                updatedCells.add( cell );
            }

            // upload update
            myWorksheet.batchUpdate( updatedCells );
        }

        return updatedServices;
    }

    private String today() {return DATE_FORMAT.format( new Date() );}

    private String todayDirectory() {return DATE_AS_DIR.format( new Date() );}

    public List<String> updateInteractionWorksheet( final SpreadsheetService service,
                                                    SpreadsheetEntry spreadsheet,
                                                    Map<String, Long> db2interactionCount ) throws Exception {
        return updateSheet( service, spreadsheet, "Interactions", db2interactionCount );
    }

    private List<String> updatePublicationWorksheet( SpreadsheetService service,
                                                     SpreadsheetEntry spreadsheet,
                                                     Map<String, Long> db2publicationsCount ) throws IOException,
            ServiceException {
        return updateSheet( service, spreadsheet, "Publications", db2publicationsCount );
    }

    //////////////////////////////////////
    // PSICQUIC/Registry related methods

    public List<PsicquicService> collectPsicquicServiceNames() throws IOException {
        /* Get list of services with metadata information from the PSICQUIC registry */
        PsicquicRegistryClient registryClient = new DefaultPsicquicRegistryClient();
        List<ServiceType> registryServices = null;
        try {
            registryServices = registryClient.listServices();
        } catch (PsicquicRegistryClientException e) {
            e.printStackTrace();
        }



        final String registryUrl = config.getPsicquicRegistryUrl();

        if ( log.isInfoEnabled() ) log.info( "Reading PSICQUIC services list from: " + registryUrl );

        // collect services via REST to to be able to pass filters on the URL (status, tags...)
        URL url = new URL( registryUrl );
        Properties props = new Properties( );
        final InputStream is = url.openStream();
        props.load( is );
        is.close();

        List<PsicquicService> services  = Lists.newArrayList();

        for ( Map.Entry<Object, Object> service : props.entrySet() ) {
            final PsicquicService s = new PsicquicService( ( String ) service.getKey(),
                    ( String ) service.getValue() );
            //todo: Refactor. This code is a bit redundant. Included because the way it was done before it was impossible to retrieve the REST url. Now the are two call to the registry
            for(ServiceType registryService: registryServices){
                if(registryService.getName().equalsIgnoreCase(s.getName())){
                    s.setRestUrl(registryService.getRestUrl());
                }
            }
            services.add( s );
            log.info( s );
        }

        if ( log.isInfoEnabled() ) log.info( "Found " + services.size() + " PSICQUIC services in the Registry." );

        return services;
    }

    /**
     * Update the given list of PSICQUIC services. that is their interaction count.
     * @param psicquicServices the list of services to update.
     */
    public Map<String, Long> updatePsicquicInteractionsStats(List<PsicquicService> psicquicServices) {
        Map<String, Long> db2interactionCount = Maps.newHashMap();

        for ( PsicquicService service : psicquicServices ) {
            log.info( service.getName() + " -> " + service.getSoapUrl() );
            PsicquicSimpleClient client = new PsicquicSimpleClient(service.getRestUrl());

            try {
                boolean error = false;

                try{
                    long count = client.countByQuery( config.getInteractionMiqlQuery() );
                    service.setInteractionCount( Long.valueOf(count).intValue() );
                    db2interactionCount.put( service.getName(), count );
                }
                catch (IOException e) {

                    int number = 1;
                    error = true;

                    while (number < numberOfTests && error == true){
                        number ++;
                        System.out.println("Failed to connect to service, try number " + number);

                        try {
                            // wait for 10 secondes and try again
                            Thread thisThread = Thread.currentThread();

                            thisThread.sleep(10000);

                            // try again

                            long count = client.countByQuery( config.getInteractionMiqlQuery() );
                            service.setInteractionCount( Long.valueOf(count).intValue() );
                            db2interactionCount.put( service.getName(), count );

                            error = false;
                        } catch (IOException e2) {
                            log.error( "An error occured while retrieving interactions from " + service.getName()+ ", test number " + number, e );
                            error = true;

                            if (number == numberOfTests){
                                log.error( "An error occured while querying PSICQUIC service: " + service.getName(), e2 );
                                sendEmail( "Failed to query " + service.getName(), ExceptionUtils.getFullStackTrace( e2 ) );
                                break;
                            }
                        }
                    }
                }
            }
            catch ( Throwable t ) {
                log.error( "An error occured while querying PSICQUIC service: " + service.getName(), t );
                sendEmail( "Failed to query " + service.getName(), ExceptionUtils.getFullStackTrace( t ) );
            }
        }

        if ( log.isDebugEnabled() ) {
            log.debug( "Interaction count collected:" );
            for ( Map.Entry<String, Long> entry : db2interactionCount.entrySet() ) {
                log.debug( entry.getKey() + " -> " + entry.getValue() );
            }
        }

        return db2interactionCount;
    }

    public Map<String, Long> collectPsicquicPublicationsStats(List<PsicquicService> psicquicServices) throws IOException {
        return collectPsicquicPublicationsStats( psicquicServices, null );
    }

    public Map<String, Long> collectPsicquicPublicationsStats( List<PsicquicService> psicquicServices,
                                                               List<String> psicquicWhiteList ) throws IOException {
        Map<String, Long> db2publicationsCount = Maps.newHashMap();

        for ( PsicquicService service : psicquicServices ) {
            if ( psicquicWhiteList != null) {
                if ( !psicquicWhiteList.contains( service.getName() ) ) {
                    // skip this resource.
                    if ( log.isInfoEnabled() ) log.info( "Skipping the count of pudmed for PSICQUIC sercice: " + service.getName() );
                    continue;
                }
            }


            int current = 0;


            if ( log.isInfoEnabled() ) log.info( "Querying PSICQUIC service: " + service.getRestUrl() );
            Set<String> pmids = Sets.newHashSet();
            boolean error = false;

            PsicquicSimpleClient simpleClient = new PsicquicSimpleClient(service.getRestUrl());
            final long totalInteractionCount = simpleClient.countByQuery( config.getPublicationMiqlQuery() );

            try {
                do {
                    PsimiTabReader mitabReader = new PsimiTabReader();
                    try {
                        Collection<BinaryInteraction> binaryInteractions = processCountPublications(current, pmids, simpleClient, mitabReader);
                        current += binaryInteractions.size();
                    } catch (IOException e) {

                        int number = 1;
                        error = true;

                        while (number < numberOfTests && error == true){
                            number ++;
                            System.out.println("Failed to connect to service, try number " + number);

                            try {
                                // wait for 10 secondes and try again
                                Thread thisThread = Thread.currentThread();

                                thisThread.sleep(10000);

                                // try again

                                Collection<BinaryInteraction> binaryInteractions = processCountPublications(current, pmids, simpleClient, mitabReader);
                                current += binaryInteractions.size();

                                error = false;
                            } catch (IOException e2) {
                                log.error( "An error occured while retrieving interactions from " + service.getName()+ ", test number " + number, e );
                                error = true;

                                if (number == numberOfTests){
                                    log.error( "An error occured while retrieving interactions from " + service.getName(), e );
                                    // email error
                                    sendEmail( "An error occured while collecting PMIDs from: " + service.getName(),
                                            ExceptionUtils.getFullStackTrace( e ) );
                                    break;
                                }
                            }
                        }
                    }

                    // show progress
                    if ( ( current % 1000 ) == 0 ) {
                        log.info( current );
                    }
                } while ( current < totalInteractionCount );
            } catch ( Throwable t ) {
                log.error( "An error occured while collecting PMIDs from " + service.getName(), t );
                error = true;

                // email error
                sendEmail( "An error occured while collecting PMIDs from: " + service.getName(),
                        ExceptionUtils.getFullStackTrace( t ) );
            }

            if ( log.isInfoEnabled() )
                log.info( "\n" + service.getName() + " -> " + pmids.size() + " publication(s)." +
                        ( error ? " However the PMID collection may have been interrupted." : "" ) );

            if ( !error ) {
                db2publicationsCount.put( service.getName(), ( long ) pmids.size() );
            } else {
                db2publicationsCount.put( service.getName(), -1L );
            }

            // Save pmid list for later processing.
            String user = null;
            if( senderEmail != null ) {
                System.out.println( "senderEmail='"+senderEmail+"'" );
                if( senderEmail.indexOf( "@" ) != -1 ) {
                    user = senderEmail.substring( 0, senderEmail.indexOf( "@" ));
                }
            }
            File parentDir = new File( (user == null ? "" : user + FILE_SEPARATOR ) + todayDirectory() );
            if ( !parentDir.exists() ) {
                parentDir.mkdirs();
            }
            writeListToDisk( pmids, new File( parentDir, service.getName() + ".pmids.txt" ) );

        } // psicquic services

        return db2publicationsCount;
    }

    private Collection<BinaryInteraction> processCountPublications(int current, Set<String> pmids, PsicquicSimpleClient simpleClient, PsimiTabReader mitabReader) throws IOException, PsimiTabException {
        InputStream result = simpleClient.getByQuery(config.getPublicationMiqlQuery(),"tab25", current, PSICQUIC_BATCH_SIZE);

        System.out.println(config.getPublicationMiqlQuery() + " / " + "tab25" + " / " + current + " / " + PSICQUIC_BATCH_SIZE);

        Collection<BinaryInteraction> binaryInteractions = mitabReader.read(result);
        boolean hasPubmed = false;

        for(BinaryInteraction binaryInteraction:binaryInteractions){
            for ( Object o : binaryInteraction.getPublications() ) {
                CrossReference cr = ( CrossReference ) o;
                if ( "pubmed".equalsIgnoreCase( cr.getDatabase() )
                        ||
                        "pmid".equalsIgnoreCase( cr.getDatabase() )) {

                    if (pmids.add( cr.getIdentifier() )){
                        System.out.println(cr.getIdentifier());
                    }

                    hasPubmed = true;
                }
            }

            if (!hasPubmed){
                log.error("Binary interaction without pubmed : " + binaryInteraction.toString());
            }
        }
        return binaryInteractions;
    }

    private void writeListToDisk( Set<String> list, File file ) throws IOException {
        BufferedWriter out = new BufferedWriter( new FileWriter( file ) );
        for ( String line : list ) {
            out.write( line + "\n" );
        }
        out.flush();
        out.close();
    }

    private String showServices( List<PsicquicService> psicquicServices ) {
        StringBuilder sb = new StringBuilder( 256 );
        sb.append( '[' );
        for ( Iterator<PsicquicService> iterator = psicquicServices.iterator(); iterator.hasNext(); ) {
            PsicquicService service = iterator.next();
            sb.append( service.getName() );
            if( iterator.hasNext() ) {
                sb.append( ", " );
            }
        }
        sb.append(']');
        return sb.toString();
    }

    //todo:separate spreedsheet update from data to be able to test data
    private void updateSpreadsheet( String email, String password, String spreadsheetKey) throws Exception {
        updateSpreadsheet(email, password, spreadsheetKey, null, null, null);
    }

    private void updateSpreadsheet( String email, String password, String spreadsheetKey, List<PsicquicService> psicquicServices,
                                    Map<String, Long> db2interactionCount, Map<String, Long> db2publicationsCount ) throws Exception {
        SpreadsheetService service = new SpreadsheetService( "PSICQUIC-stats-collector-1" );
        // 2010-02: necessary or batch update fail with error like: "If-Match or If-None-Match header required"
        service.setProtocolVersion( SpreadsheetService.Versions.V1 );
        service.setUserCredentials( email, password );

        final SpreadsheetEntry spreadsheetEntry = SpreadsheetFacade.getSpreadsheetWithKey( service, spreadsheetKey );
        if(psicquicServices == null){
            psicquicServices = collectPsicquicServiceNames();
        }

        // Update interaction counts
        if(db2interactionCount == null){
            db2interactionCount = updatePsicquicInteractionsStats( psicquicServices ); //data
        }
        List<String> updatedServices = updateInteractionWorksheet( service, spreadsheetEntry, db2interactionCount ); //spreadsheet
        if ( log.isInfoEnabled() ) log.info( updatedServices.size() + " services updated: " + updatedServices );

        // Update publication for those services that have a different count of interactions
        if( ! updatedServices.isEmpty()) {
            if(db2publicationsCount == null){
                db2publicationsCount = collectPsicquicPublicationsStats( psicquicServices, updatedServices); //data
            }
            updatedServices.addAll(updatePublicationWorksheet(service, spreadsheetEntry, db2publicationsCount)); //spreadsheet
        } else {
            log.info( "No services with interaction count update, skipping publication update." );
        }

        if( ! updatedServices.isEmpty() ) {
            sendEmail( "Spreadsheet was updated", "Some services had new data online: " + showServices(psicquicServices) );
        } else {
            sendEmail( "Spreadsheet was NOT updated", "No PSICQUIC services had new data." );
        }
    }

    /////////////////
    // M A I N

    public static void main( String[] args ) throws Exception {

        // TODO Auto resize of the spreadsheet when we reach the maximum size

        if ( args.length < 3 ) {
            System.err.println( "usage: PsicquicStatsCollector <gmail.account> <password> <spreadsheet.key> " +
                    "[-D"+PSICQUIC_REGISTRY_URL_KEY+"=pathToFile] " +
                    "[-D"+SMTP_CONFIG_FILE_KEY+"=pathToFile]" );
            System.exit( 1 );
        }

        // Handle input parameters
        final String email = args[0];
        final String password = args[1];
        final String spreadsheetKey = args[2];

        PsicquicStatsCollector collector = new PsicquicStatsCollector();
        collector.updateSpreadsheet( email, password, spreadsheetKey );
    }




}
