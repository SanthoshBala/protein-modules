package org.hupo.psi.mi.psicquic.clustering;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.stereotype.Component;

import java.io.File;
import java.io.Serializable;

/**
 * Clustering config.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 1.0
 */
@Component
public class ClusteringConfig implements Serializable {

    private static final Log log = LogFactory.getLog( ClusteringConfig.class );

    public static final String SYSTEM_KEY_DATA_LOCATION = "psicquic.clustering.data";

    public static final String TEMP_DIR = System.getProperty( "java.io.tmpdir" );
    public static final String FILE_SEPARATOR = System.getProperty( "file.separator" );

    /**
     * Fully qualified filename of the directory that stores the computed data.
     */
    private String dataLocation;

    public ClusteringConfig() {
    }

    private String buildDefaultDir() {
        String defaultLocation = TEMP_DIR;
        if( ! defaultLocation.endsWith( FILE_SEPARATOR ) ) {
            defaultLocation += FILE_SEPARATOR;
        }
        defaultLocation += "clustering-data";
        createDirectoryIfNotExist( defaultLocation );
        return defaultLocation;
    }

    private File createDirectoryIfNotExist( String dir ) {
        File file = new File( dir );
        if(!file.exists()) {
           file.mkdirs();
        }
        return file;
    }

    public String getDataLocation() {
        if( dataLocation == null ) {
            dataLocation = buildDefaultDir();
        }
        return dataLocation;
    }

    public File getDataLocationFile() {
        return new File(getDataLocation());
    }

    public void setDataLocation( String dataLocation ) {
        File file = createDirectoryIfNotExist( dataLocation );
        this.dataLocation = file.getAbsolutePath();
        log.warn( "Setting data location to: " + this.dataLocation );
    }
}
