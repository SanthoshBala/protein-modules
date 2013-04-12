package org.hupo.psi.mi.psicquic.clustering.job;

import org.hupo.psi.mi.psicquic.clustering.Service;

import java.io.UnsupportedEncodingException;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 * Generates random string for our ClusteringJob ids.
 * <p/>
 * Logic of this generator based on: http://stackoverflow.com/questions/41107/how-to-generate-a-random-alpha-numeric-string-in-java
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class JobIdGenerator {

    private SecureRandom random = new SecureRandom();

    public JobIdGenerator() {
    }

    public String getNextJobId() {
        // random id
        return new BigInteger(130, random).toString(32);
    }

    private String convertToHex(byte[] data) {
        StringBuffer buf = new StringBuffer();
        for (int i = 0; i < data.length; i++) {
            int halfbyte = (data[i] >>> 4) & 0x0F;
            int two_halfs = 0;
            do {
                if ((0 <= halfbyte) && (halfbyte <= 9))
                    buf.append((char) ('0' + halfbyte));
                else
                    buf.append((char) ('a' + (halfbyte - 10)));
                halfbyte = data[i] & 0x0F;
            } while(two_halfs++ < 1);
        }
        return buf.toString();
    }

    private String sha1(String text) throws NoSuchAlgorithmException, UnsupportedEncodingException {
        MessageDigest md;
        md = MessageDigest.getInstance("SHA-1");
        byte[] sha1hash = new byte[40];
        md.update(text.getBytes("iso-8859-1"), 0, text.length());
        sha1hash = md.digest();
        return convertToHex(sha1hash);
    }

    public String generateJobId( ClusteringJob job ) {
        StringBuilder sb = new StringBuilder(512);
        sb.append( job.getMiql().trim() ).append( "/" );

        // sort the list of services alphabetically
        for ( Service service : job.getServices() ) {
            sb.append( service.getName() ).append( "_" );
        }

        // hash the result
        try {
            return sha1( sb.toString() );
        } catch ( Exception e ) {
            throw new IllegalStateException( "Could not generated a jobId based on job: " + job, e );
        }
    }
}
