package playground.dbpairs;

import psidev.psi.mi.tab.model.BinaryInteraction;

import java.io.Serializable;
import java.util.*;

public class Report implements Serializable {

    // TODO currently the counts are representing instances (potentially multiple per MITAB lines),
    // TODO it would be nice to be able to count MITAB lines instead

    // databases collected in column identifier of molecule A and B
    private Set<String> databases4idA = new HashSet<String>();
    private Set<String> databases4idB = new HashSet<String>();

    // databases collected in column alternative identifier of molecule A and B
    private Set<String> databases4AltidA = new HashSet<String>();
    private Set<String> databases4AltidB = new HashSet<String>();

    private Set<String> aliases4A = new HashSet<String>();
    private Set<String> aliases4B = new HashSet<String>();
    private Set<String> aliases = new HashSet<String>();

    // databases collected in both columns identifier and alternative identifier of molecule A and B
    private Set<String> databases4id = new HashSet<String>();

    private Collection<String> databasePairs = new HashSet<String>();


    // Instance Counter

    // count usage of a specific database in column identifier of molecule A and B
    private Map<String, Integer> idDatabase2count = new HashMap<String, Integer>();

    // count usage of a specific database in column alternative identifier of molecule A and B
    private Map<String, Integer> altIdDatabase2count = new HashMap<String, Integer>();

    // count usage of a specific database in both columns identifier and alternative identifier of molecule A and B
    private Map<String, Integer> database2count = new HashMap<String, Integer>();

    // count of usage of specific database/type from the alias column of interactor A and B
    private Map<String, Integer> alias2count = new HashMap<String, Integer>();

    private Map<String, Integer> databasePairs2count = new HashMap<String, Integer>();


    // MITAB Line counter

    // these sets are used to make sure that we only count each instance once per MITAB line
    private Set<String> id4line = new HashSet<String>();
    private Set<String> altId4line = new HashSet<String>();
    private Set<String> aliases4line = new HashSet<String>();
    private Set<String> dbPairs4line = new HashSet<String>();

    private Map<String, Integer> idDatabase2lineCount = new HashMap<String, Integer>();
    private Map<String, Integer> altIdDatabase2lineCount = new HashMap<String, Integer>();
    private Map<String, Integer> database2lineCount = new HashMap<String, Integer>();
    private Map<String, Integer> alias2lineCount = new HashMap<String, Integer>();
    private Map<String, Integer> databasePairs2lineCount = new HashMap<String, Integer>();


    private int lineCount = 0;


    private void countInstance( Map<String, Integer> db2count, String db ) {
        if( !db2count.containsKey( db ) ) {
            db2count.put( db, 1 );
        } else {
            int count = db2count.get( db );
            db2count.put( db, count + 1 );
        }
    }


    //////////////////
    // Report feeding

    public void addDatabaseA( String db ) {
        databases4idA.add( db );
        databases4id.add( db );
        countInstance( idDatabase2count, db );
        countInstance( database2count, db );

        // done only once for a given db in the current MITAB line
        if( ! id4line.contains( db ) ) {
            id4line.add( db );
            countInstance( idDatabase2lineCount, db );
            countInstance( database2lineCount, db );
        }
    }

    public void addDatabaseB( String db ) {
        databases4idB.add( db );
        databases4id.add( db );
        countInstance( idDatabase2count, db );
        countInstance( database2count, db );

        // done only once for a given db in the current MITAB line
        if( ! id4line.contains( db ) ) {
            id4line.add( db );
            countInstance( idDatabase2lineCount, db );
            countInstance( database2lineCount, db );
        }
    }

    public void addAlternativeDatabaseA( String db ) {
        databases4AltidA.add( db );
        databases4id.add( db );
        countInstance( altIdDatabase2count, db );
        countInstance( database2count, db );

        // done only once for a given db in the current MITAB line
        if( ! altId4line.contains( db ) ) {
            altId4line.add( db );
            countInstance( altIdDatabase2lineCount, db );
            countInstance( database2lineCount, db );
        }
    }

    public void addAlternativeDatabaseB( String db ) {
        databases4AltidB.add( db );
        databases4id.add( db );
        countInstance( altIdDatabase2count, db );
        countInstance( database2count, db );

        // done only once for a given db in the current MITAB line
        if( ! altId4line.contains( db ) ) {
            altId4line.add( db );
            countInstance( altIdDatabase2lineCount, db );
            countInstance( database2lineCount, db );
        }
    }

    public void addAliasA( String alias ) {
        aliases4A.add( alias );
        aliases.add( alias );
        countInstance( alias2count, alias );

        // done only once for a given db in the current MITAB line
        if( ! aliases4line.contains( alias ) ) {
            aliases4line.add( alias );
            countInstance( alias2lineCount, alias );
        }
    }

    public void addAliasB( String alias ) {
        aliases4B.add( alias );
        aliases.add( alias );
        countInstance( alias2count, alias );

        // done only once for a given db in the current MITAB line
        if( ! aliases4line.contains( alias ) ) {
            aliases4line.add( alias );
            countInstance( alias2lineCount, alias );
        }
    }

    public void addDatabasePairs( Collection<String> dbA, Collection<String> dbB ) {
        for ( String a : dbA ) {
            for ( String b : dbB ) {
                String pair;
                if( a.compareTo( b ) < 0 ) {
                   pair = a+"-"+b;
                } else {
                    pair = b+"-"+a;
                }
                databasePairs.add( pair );
                countInstance( databasePairs2count, pair );

                // done only once for a given db in the current MITAB line
                if( ! dbPairs4line.contains( pair ) ) {
                    dbPairs4line.add( pair );
                    countInstance( databasePairs2lineCount, pair );
                }
            }
        }
    }

    public void incrementLineCount(){
        lineCount++;
    }


    ///////////////////////////
    // Getters and Setters

    public Collection<String> getAllDatabases() {
        return databases4id;
    }

    public Collection<String> getDatabasePairs() {
        return databasePairs;
    }

    public Map<String, Integer> getIdDatabase2count() {
        return idDatabase2count;
    }

    public Map<String, Integer> getAltIdDatabase2count() {
        return altIdDatabase2count;
    }

    public Map<String, Integer> getDatabase2count() {
        return database2count;
    }

    public Map<String, Integer> getDatabasePairs2count() {
        return databasePairs2count;
    }

    public Map<String, Integer> getAlias2count() {
        return alias2count;
    }

    public Map<String, Integer> getIdDatabase2lineCount() {
        return idDatabase2lineCount;
    }

    public Map<String, Integer> getAltIdDatabase2lineCount() {
        return altIdDatabase2lineCount;
    }

    public Map<String, Integer> getDatabase2lineCount() {
        return database2lineCount;
    }

    public Map<String, Integer> getAlias2lineCount() {
        return alias2lineCount;
    }

    public Map<String, Integer> getDatabasePairs2lineCount() {
        return databasePairs2lineCount;
    }

    public int getLineCount() {
        return lineCount;
    }


    //////////////////////////////
    // MITAB line scope control

    public void startMitabLine( BinaryInteraction bi ) {
        id4line.clear();
        altId4line.clear();
        aliases4line.clear();
        dbPairs4line.clear();
    }

    public void endMitabLine( BinaryInteraction bi ) {

    }
}
