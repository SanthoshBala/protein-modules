package playground.dbpairs;

import org.apache.commons.lang.StringUtils;
import playground.dbpairs.selectors.*;
import psidev.psi.mi.tab.PsimiTabReader;
import psidev.psi.mi.tab.model.Alias;
import psidev.psi.mi.tab.model.BinaryInteraction;
import psidev.psi.mi.tab.model.CrossReference;

import java.io.*;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class AnalyseMitab {

    public static final DecimalFormat NUMBER_FORMATTER = new DecimalFormat( "###,###" );

    public static final DecimalFormat PERCENT_FORMATTER = new DecimalFormat( "###.##" );

    public static final SimpleDateFormat DATE_FORMATTER = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss");


    ///////////////////////////
    // Private methods

    private Map<String, Collection<String>> buildItem2Providers( final Map<String, Report> provider2report,
                                                                 final CollectionSelector selector ) {

        final Map<String, Collection<String>> item2providers = new HashMap<String, Collection<String>>();

        for ( Map.Entry<String, Report> entry : provider2report.entrySet() ) {
            final String provider = entry.getKey();
            final Report report = entry.getValue();

            final Collection<String> collection = selector.select( report );
            for ( String item : collection ) {
                if ( !item2providers.containsKey( item ) ) {
                    item2providers.put( item, new ArrayList<String>() );
                }
                final Collection<String> providers = item2providers.get( item );
                providers.add( provider );
            }
        }

        return item2providers;
    }

    private Map<String, Collection<String>> buildItem2ProviderCount( final Map<String, Report> provider2report,
                                                                     final MapSelector selector,
                                                                     final boolean printLinePercentage) {

        final Map<String, Collection<String>> item2providers = new HashMap<String, Collection<String>>();

        for ( Map.Entry<String, Report> entry : provider2report.entrySet() ) {
            final String provider = entry.getKey();
            final Report report = entry.getValue();

            final Map<String, Integer> items2count = selector.select( report );
            for ( Map.Entry<String, Integer> item2count : items2count.entrySet() ) {
                final String item = item2count.getKey();
                final Integer count = item2count.getValue();

                if ( !item2providers.containsKey( item ) ) {
                    item2providers.put( item, new ArrayList<String>() );
                }
                final Collection<String> providers = item2providers.get( item );
//                providers.add( provider + "(" + NUMBER_FORMATTER.format(count) + ")" );

                String percentStr = null;
                if( printLinePercentage ) {
                    final int lineCount = report.getLineCount();
                    float percent = ((float)count / lineCount ) * 100;
                    percentStr = PERCENT_FORMATTER.format( percent ) + "%";
                }
//                providers.add( provider + "(" + NUMBER_FORMATTER.format(count) + (percentStr != null ? percentStr : "") + ")" );
                providers.add( NUMBER_FORMATTER.format(count) + " lines, i.e. " + (percentStr != null ? percentStr : "") );
            }
        }

        return item2providers;
    }

    private Map<String, Collection<String>> buildProvider2ItemCount( Map<String, Report> provider2report,
                                                                     MapSelector selector,
                                                                     boolean printLinePercentage ) {

        final Map<String, Collection<String>> provider2items = new HashMap<String, Collection<String>>();

        for ( Map.Entry<String, Report> entry : provider2report.entrySet() ) {
            final String provider = entry.getKey();
            final Report report = entry.getValue();

            final Map<String, Integer> items2count = selector.select( report );
            for ( Map.Entry<String, Integer> item2count : items2count.entrySet() ) {
                final String item = item2count.getKey();
                final Integer count = item2count.getValue();

                if ( !provider2items.containsKey( provider ) ) {
                    provider2items.put( provider, new ArrayList<String>() );
                }
                final Collection<String> items = provider2items.get( provider );
                String percentStr = null;
                if( printLinePercentage ) {
                    final int lineCount = report.getLineCount();
                    float percent = ((float)count / lineCount ) * 100;
                    percentStr = ", " + PERCENT_FORMATTER.format( percent ) + "%";
                }
                items.add( item + "(" + NUMBER_FORMATTER.format(count) + (percentStr != null ? percentStr : "") + ")" );
            }
        }

        return provider2items;
    }

    private Report analyse( File file ) throws Exception {
        return analyse( file, -1 ); // no limit
    }

    private Report analyse( File file, final int maxLines ) throws Exception {

        boolean hasLineCountLimit = ( maxLines != -1 );

        final Report report = new Report();
        final PsimiTabReader reader = new PsimiTabReader( false );
        final Iterator<BinaryInteraction> iterator = reader.iterate( file );

        int lineCount = 0;
        while ( iterator.hasNext() ) {

            if( hasLineCountLimit ) {
                if( lineCount >= maxLines ) {
                    break;
                }
            }

            lineCount++;

            final BinaryInteraction bi = iterator.next();
            report.startMitabLine( bi );
            report.incrementLineCount();

            Collection<String> dbIdentifierA = getDatabases( bi.getInteractorA().getIdentifiers() );
            for ( String database : dbIdentifierA ) {
                report.addDatabaseA( database );
            }

            Collection<String> dbIdentifierB = getDatabases( bi.getInteractorB().getIdentifiers() );
            for ( String database : dbIdentifierB ) {
                report.addDatabaseB( database );
            }

            Collection<String> dbAltIdentifierA = getDatabases( bi.getInteractorA().getAlternativeIdentifiers() );
            for ( String database : dbAltIdentifierA ) {
                report.addAlternativeDatabaseA( database );
            }

            Collection<String> dbAltIdentifierB = getDatabases( bi.getInteractorB().getAlternativeIdentifiers() );
            for ( String database : dbAltIdentifierB ) {
                report.addAlternativeDatabaseB( database );
            }

            Collection<String> aliasesA = getAliases( bi.getInteractorA().getAliases() );
            for ( String a : aliasesA ) {
                report.addAliasA( a );
            }

            Collection<String> aliasesB = getAliases( bi.getInteractorB().getAliases() );
            for ( String a : aliasesB ) {
                report.addAliasB( a );
            }

            // Build DB pairs, based on both identifiers and alternative identifier between molecule A and B
            final Set<String> allIdsA = new HashSet<String>();
            allIdsA.addAll( dbIdentifierA );
            allIdsA.addAll( dbAltIdentifierA );

            final Set<String> allIdsB = new HashSet<String>();
            allIdsB.addAll( dbIdentifierB );
            allIdsB.addAll( dbAltIdentifierB );

            report.addDatabasePairs( allIdsA, allIdsB );

            report.endMitabLine( bi );
        }

        System.out.println( "  > " + NUMBER_FORMATTER.format(lineCount) + " MITAB lines processed.");

        return report;
    }

    private Collection<String> getAliases( Collection<Alias> aliases ) {
        Set<String> uniqAliases = new HashSet<String>();
        for ( Alias alias : aliases ) {
            StringBuilder sb = new StringBuilder();
            sb.append( alias.getDbSource() );
            if ( alias.getAliasType() != null ) {
                sb.append( '/' ).append( alias.getAliasType() );
            }
            uniqAliases.add( sb.toString() );
        }
        return uniqAliases;
    }

    private Collection<String> getDatabases( Collection<CrossReference> crossReferences ) {
        Set<String> dbs = new HashSet<String>();
        for ( CrossReference cr : crossReferences ) {
            dbs.add( cr.getDatabase() );
        }
        return dbs;
    }

    private void printThings2Providers( Map<String, Collection<String>> thing2provider,
                                        boolean showSingleton,
                                        String fromHeader,
                                        String toHeader ) {

        if( ! hasData( thing2provider, showSingleton ) ) {
            System.out.println( "No data." );
            return;
        }

        if( ! showSingleton ) {
            System.out.println( "WARNING - The table below is only showing data involving more than one provider.\n" );
        }

        int longestKey = findLongestKey( thing2provider );

        // header
        System.out.print( StringUtils.rightPad( fromHeader, longestKey + 2 ) );
        System.out.println( toHeader );

        // header separator
        System.out.print( StringUtils.rightPad( "", longestKey, '-' ) );
        System.out.print( "  " );
        System.out.println( StringUtils.rightPad( "", 40, '-' ) );

        // data
        for ( Map.Entry<String, Collection<String>> entry : thing2provider.entrySet() ) {
            final String thing = entry.getKey();
            final Collection<String> providers = entry.getValue();
            if ( providers.size() > 1 ) {
                System.out.println( StringUtils.rightPad( thing, longestKey + 2 ) + providers );
            } else {
                if( showSingleton ) {
                    System.out.println( StringUtils.rightPad( thing, longestKey + 2 ) + providers );
                }
            }
        }
    }

    private boolean hasData( Map<String, Collection<String>> thing2provider, boolean showSingleton ) {
        final int minInstance = (showSingleton ? 1 : 2);
        for ( Map.Entry<String, Collection<String>> entry : thing2provider.entrySet() ) {
            if( entry.getValue().size() >= minInstance ) {
                return true;
            }
        }
        return false;
    }

    private int findLongestKey( final Map<String, Collection<String>> map ) {
        int max = 0;
        for ( String key : map.keySet() ) {
            max = Math.max( max, key.length() );
        }
        return max;
    }

    private void printHeader( final String headerTitle, String s ) {
        final int headerWidth = headerTitle.length() + 20;
        System.out.println("\n\n");
        System.out.println( StringUtils.rightPad( "", headerWidth, s ) );
        System.out.println( StringUtils.center( headerTitle, headerWidth, s ) );
        System.out.println();
    }


    /////////////////
    // Public methods

    public void runAnalysis( final CountingStrategy countingStrategy ) {
        runAnalysis( -1, countingStrategy );
    }

    public void runAnalysis( final int maxLines, final CountingStrategy countingStrategy ) {

        Map<String, Report> provider2report = new HashMap<String, Report>();

//        File dir = new File( "/Users/skerrien-old/psicquic-imex-data" );
//        final String[] files = dir.list( new FilenameFilter() {
//            public boolean accept( File dir, String name ) {
//                return name.endsWith( ".tsv" );
//            }
//        } );
//
//        System.out.println( "Found " + files.length + " files." );
//
//        if( maxLines != -1 ) {
//            System.out.println("\nWARNING - The current analysis is only reading the first " + maxLines + " MITAB lines of each provider.\n");
//        }
//
//        for ( int i = 0; i < files.length; i++ ) {
//            String file = files[i];
//            System.out.println( "Processing " + file );
//
//            try {
//                Report report = analyse( new File( dir, file ), maxLines );
//                String providerName = file.replaceAll( "\\.tsv", "" );
//                provider2report.put( providerName, report );
//            } catch ( Exception e ) {
//                System.err.println( "Count not parse: " + file );
//                e.printStackTrace();
//            }
//        } // all files
//
//
//        /////////////////////
//        // Save current report
//
//        if( ! provider2report.isEmpty() ) {
//            String filename = DATE_FORMATTER.format( new Date() ) + ".ser";
//            try {
//                saveReports( provider2report, new File( filename ) );
//            } catch ( IOException e ) {
//                System.err.println( "Failed to save results to disk." );
//                e.printStackTrace();
//            }
//        }

        try {
            provider2report = loadReports( new File( "2011.03.30.13.16.26.ser" ) );
        } catch ( Exception e ) {
            e.printStackTrace();
        }

        /////////////////////////
        // Print Results

//        printAnalysisResults( provider2report, countingStrategy );
        printAnalysisResultByProviders( provider2report, countingStrategy );
    }

    public void printAnalysisResults( Map<String, Report> provider2report, CountingStrategy countingStrategy ) {
        printHeader( " Usage of database (identifier column) ", "=" );
        final boolean showLinePercentage = countingStrategy.isShowingLinePercentage();

        final MapSelector idDatabaseSelector = countingStrategy.getIdDatabaseSelector();
        Map<String, Collection<String>> idDatabase2providerCount = buildItem2ProviderCount( provider2report,
                                                                                            idDatabaseSelector,
                                                                                               showLinePercentage );
        printThings2Providers( idDatabase2providerCount, true, "Identifier", "Providers" );

        if( provider2report.size() > 1 ) {
            Map<String, Collection<String>> provider2IdDatabasesCount = buildProvider2ItemCount( provider2report,
                                                                                                 idDatabaseSelector,
                                                                                                 showLinePercentage);

            System.out.println( "\n----  Breakdown per provider  ----\n" );
            printThings2Providers( provider2IdDatabasesCount, true, "Provider", "Identifiers" );
        }


        printHeader( " Usage of database (alt. identifier column) ", "=" );
        final MapSelector altIdDatabaseSelector = countingStrategy.getAltIdDatabaseSelector();
        Map<String, Collection<String>> altIdDatabase2providerCount = buildItem2ProviderCount( provider2report,
                                                                                               altIdDatabaseSelector,
                                                                                               showLinePercentage );
        printThings2Providers( altIdDatabase2providerCount, true, "Identifier", "Providers" );

        if( provider2report.size() > 1 ) {
            System.out.println( "\n----  Breakdown per provider  ----\n" );
            Map<String, Collection<String>> provider2altIdDatabasesCount = buildProvider2ItemCount( provider2report,
                                                                                                    altIdDatabaseSelector,
                                                                                                    showLinePercentage );
            printThings2Providers( provider2altIdDatabasesCount, true, "Provider", "Identifiers" );
        }


        printHeader( " Usage of aliases ", "=" );
        final MapSelector aliasesSelector = countingStrategy.getAliasSelector();
        Map<String, Collection<String>> alias2providerCount = buildItem2ProviderCount( provider2report,
                                                                                       aliasesSelector,
                                                                                               showLinePercentage );
        printThings2Providers( alias2providerCount, true, "Alias(db/type)", "Providers" );

        if( provider2report.size() > 1 ) {
            System.out.println( "\n----  Breakdown per provider  ----\n" );
            Map<String, Collection<String>> provider2aliasesCount = buildProvider2ItemCount( provider2report,
                                                                                             aliasesSelector,
                                                                                             showLinePercentage );
            printThings2Providers( provider2aliasesCount, true, "Provider", "Aliases(db/type)" );
        }


        printHeader( " Usage of database-pairs ", "=" );
        final MapSelector pairsSelector = countingStrategy.getDatabasePairsSelector();
        Map<String, Collection<String>> databasePair2providersCount = buildItem2ProviderCount( provider2report,
                                                                                               pairsSelector,
                                                                                               showLinePercentage );
        printThings2Providers( databasePair2providersCount, true, "Database Pair", "Providers" );

        if( provider2report.size() > 1 ) {
            System.out.println( "\n\n----  Breakdown per provider  ----\n" );
            Map<String, Collection<String>> provider2databasePairsCount = buildProvider2ItemCount( provider2report,
                                                                                                   pairsSelector,
                                                                                                   showLinePercentage );
            printThings2Providers( provider2databasePairsCount, true, "Provider", "Database Pairs" );
        }
    }

    public void printAnalysisResultByProviders( Map<String, Report> provider2report, CountingStrategy countingStrategy ) {

        for ( final String provider : provider2report.keySet() ) {
            final Report report = provider2report.get( provider );

            printHeader( provider, "#" );

            final Map<String, Report> singleProvider2report = new HashMap<String, Report>();
            singleProvider2report.put( provider,  report );

            printAnalysisResults(singleProvider2report, countingStrategy );
        }
    }

    public void printAnalysisResults( File backupFile, CountingStrategy countingStrategy ) throws ClassNotFoundException, IOException {
        printAnalysisResults( loadReports( backupFile ), countingStrategy );
    }

    public void saveReports( Map<String, Report> provider2report, File backupFile ) throws IOException {
        System.out.println("\nSaving reports in " + backupFile.getAbsolutePath() );
        ObjectOutput out = new ObjectOutputStream(new FileOutputStream( backupFile ));
        out.writeObject(provider2report);
        out.close();
        System.out.println("done.");
    }

    public Map<String, Report> loadReports( File backupFile ) throws IOException, ClassNotFoundException {
        System.out.println("\nLoading reports from " + backupFile.getAbsolutePath() );
        ObjectInputStream in = new ObjectInputStream(new FileInputStream(backupFile));
        Map<String, Report> provider2report = (Map<String, Report>) in.readObject();
        in.close();
        return provider2report;
    }

    //////////////////
    // DEMO

    public static void main( String[] args ) {


        // TODO add the type of the xref when available in the stats

        final AnalyseMitab analyseMitab = new AnalyseMitab();
        analyseMitab.runAnalysis( new LineCountingStrategy() );

    }
}
