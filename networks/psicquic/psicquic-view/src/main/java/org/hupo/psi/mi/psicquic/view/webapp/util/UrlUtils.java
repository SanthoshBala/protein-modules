package org.hupo.psi.mi.psicquic.view.webapp.util;

/**
 * URL related utilities.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 1.2
 */
public class UrlUtils {

    private UrlUtils() {
    }

    /**
     * Encode non ASCII characters in URLs into ASCII.
     * <p/>
     * Reference: http://www.w3schools.com/TAGS/ref_urlencode.asp
     * @param url
     * @return
     */
    public static String encodeUrl(String url) {
        if( url == null ) return url;

        url = url.replaceAll( " ", "%20" );
        url = url.replaceAll( "\"", "%22" );
        url = url.replaceAll( "\\(", "%28" );
        url = url.replaceAll( "\\)", "%29" );
        url = url.replaceAll( ":", "%3A" );

        return url;
    }
}
