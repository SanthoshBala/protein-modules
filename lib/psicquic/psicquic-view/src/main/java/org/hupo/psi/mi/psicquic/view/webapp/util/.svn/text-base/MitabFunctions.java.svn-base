/**
 * Copyright 2008 The European Bioinformatics Institute, and others.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hupo.psi.mi.psicquic.view.webapp.util;

import org.apache.commons.lang.StringUtils;
import psidev.psi.mi.tab.converter.tab2graphml.AliasByTypeComparator;
import psidev.psi.mi.tab.converter.tab2graphml.AliasComparator;
import psidev.psi.mi.tab.model.Alias;
import psidev.psi.mi.tab.model.CrossReference;
import psidev.psi.mi.tab.model.Interactor;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.util.*;

/**
 * Functions to be used in the UI to control the display.
 *
 * @author Prem Anand (prem@ebi.ac.uk)
 * @version $Id$
 * @since 2.0.1
 */
public final class MitabFunctions {

    private static final String ENZYME_PSI_REF = "MI:0501";

    private static final String ENZYME_TARGET_PSI_REF = "MI:0502";

    private static final String DRUG_PSI_REF = "MI:1094";

    private static final String DRUG_TARGET_PSI_REF = "MI:1095";


    private MitabFunctions() {
    }

    public static String getInteractorDisplayName(Interactor interactor) {
         return interactor.getIdentifiers().iterator().next().getIdentifier();
    }

    public static boolean isEnzyme( Collection<CrossReference> crossReferences ) {

        for ( CrossReference crossReference : crossReferences ) {
            if ( ENZYME_PSI_REF.equals( crossReference.getIdentifier() ) ) {
                return true;
            }
        }
        return false;
    }

    public static boolean isEnzymeTarget( Collection<CrossReference> crossReferences ) {

        for ( CrossReference crossReference : crossReferences ) {
            if ( ENZYME_TARGET_PSI_REF.equals( crossReference.getIdentifier() ) ) {
                return true;
            }
        }
        return false;
    }


    public static boolean isDrug( Collection<CrossReference> crossReferences ) {

        for ( CrossReference crossReference : crossReferences ) {
            if ( DRUG_PSI_REF.equals( crossReference.getIdentifier() ) ) {
                return true;
            }
        }
        return false;
    }

    public static boolean isDrugTarget( Collection<CrossReference> crossReferences ) {

        for ( CrossReference crossReference : crossReferences ) {
            if ( DRUG_TARGET_PSI_REF.equals( crossReference.getIdentifier() ) ) {
                return true;
            }
        }
        return false;
    }



    public static String[] getInitialForCrossReference( Collection<CrossReference> crossReferences ) {

        String[] roleAndDesc = new String[2];
        String role;
        if ( crossReferences.size() == 1 ) {
            role = crossReferences.iterator().next().getText();
            if ( role != null && !( "unspecified role".equals( role ) ) ) {
                roleAndDesc[0] = getFirstLetterofEachToken( role ).toUpperCase();
                roleAndDesc[1] = roleAndDesc[0] + " = " + role;
            } else {
                roleAndDesc[0] = "-";
                roleAndDesc[1] = "-";
            }
        } else {
            Set<String> rolesSymbol = new HashSet<String>();
            Set<String> rolesDesc = new HashSet<String>();

            for ( CrossReference crossReference : crossReferences ) {
                if ( crossReference.getText() != null && !"unspecified role".equals( crossReference.getText() ) ) {
                    String symbol = getFirstLetterofEachToken( crossReference.getText());
                    rolesSymbol.add( symbol );
                    rolesDesc.add( symbol + " = "+ crossReference.getText() );
                }
            }

            if ( rolesSymbol.size() == 1 ) {
                String symbol = rolesSymbol.iterator().next();
                String desc = rolesDesc.iterator().next();
                if ( symbol != null && symbol.length() > 0 ) {
                    roleAndDesc[0] = symbol;
                    roleAndDesc[1] = symbol + " = "+desc;
                } else {
                    roleAndDesc[0] = "-";
                    roleAndDesc[1] = "-";
                }
            }
            else if ( rolesSymbol.size() > 1 ) {
                roleAndDesc[0] = StringUtils.join( rolesSymbol.toArray(), "," );
                roleAndDesc[1] = StringUtils.join( rolesDesc.toArray(), "," );
            }

        }
        return roleAndDesc;
    }


    public static String getFirstLetterofEachToken(String stringToken) {
        String s = "";
        if (stringToken.split("\\s+").length == 1) {
            return stringToken.substring(0, 2).toUpperCase();
        }
        for (String str : stringToken.split("\\s+")) {
            s = s + str.substring(0, 1).toUpperCase();
        }
        return s;
    }

    public static String getIntactIdentifierFromCrossReferences(Collection<CrossReference> xrefs) {
        return getIdentifierFromCrossReferences(xrefs, "intact");
    }

     public static String getUniprotIdentifierFromCrossReferences(Collection<CrossReference> xrefs) {
        return getIdentifierFromCrossReferences(xrefs, "uniprotkb");
    }

    public static String getChebiIdentifierFromCrossReferences(Collection<CrossReference> xrefs) {
        return getIdentifierFromCrossReferences(xrefs, "chebi");
    }

    public static String getIdentifierFromCrossReferences(Collection<CrossReference> xrefs, String databaseLabel) {
        for (CrossReference xref : xrefs) {
            if (databaseLabel.equals(xref.getDatabase())) {
                return xref.getIdentifier();
            }
        }
        return null;
    }

    public static boolean isApprovedDrug( String drugType ) {
        if ( drugType != null ) {
            if ( drugType.toLowerCase().contains( "approved".toLowerCase() ) ) {
                return true;
            }
        }
        return false;
    }

    public static Collection getFilteredCrossReferences( Collection<CrossReference> xrefs, String filter ) {
        if ( filter == null ) {
//            throw new NullPointerException( "You must give a non null filter" );
			return xrefs;
        }

        List<CrossReference> filteredList = new ArrayList<CrossReference>();

        for ( CrossReference xref : xrefs ) {
            if ( filter.equals( xref.getText() ) ) {
                filteredList.add( xref );
            }
        }
        return filteredList;
    }

	public static Collection getFilteredCrossReferencesByChecksum( Collection<CrossReference> xrefs) {
		List<CrossReference> filteredList = new ArrayList<CrossReference>();

//		final List<String> checksumNames = new ArrayList<String>(Arrays.asList(new String[]
//				{"checksum", "smiles string", "standard inchi", "inchi key",
//						"standard inchi key", "rogid", "rigid", "crogid", "crc", "crc64"}));

		final List<String> checksumNames = new ArrayList<String>(Arrays.asList(new String[]
				{"checksum", "rogid", "rigid", "crogid", "crc", "crc64"}));

		final String IREFINDEX = "irefindex";


		for ( CrossReference xref : xrefs ) {
			String aux = xref.getDatabase();
			if (!checksumNames.contains(aux) && !IREFINDEX.equalsIgnoreCase(aux)) {
				filteredList.add( xref );
			}
		}
		return filteredList;
	}

    /**
     * Filter the given collection by removing any xref that have any of the two filters.
     * @param xrefs
     * @param databaseFilter
     * @param textFilter
     * @return
     */
    public static Collection getExclusionFilteredCrossReferences( Collection<CrossReference> xrefs, String textFilter, String databaseFilter) {
        if ( databaseFilter == null && textFilter == null) {
//            throw new NullPointerException( "You must give at least one non null filter" );
			return xrefs;
        }


        List<CrossReference> filteredList = new ArrayList<CrossReference>();

        for ( CrossReference xref : xrefs ) {
            if ((databaseFilter!=null && !databaseFilter.equals( xref.getDatabase())) && (textFilter!=null && !textFilter.equals(xref.getText()))){
                filteredList.add( xref );
            }
        }
        return filteredList;
    }

	public static Collection getExclusionFilteredAliases( Collection<Alias> aliases ) {

		List<Alias> filteredList = new ArrayList<Alias>();

		for ( Alias alias : aliases ) {
			if ( alias.getAliasType() == null ||
					(alias.getAliasType() != null && !alias.getAliasType().equals("display_short") && !alias.getAliasType().equals("display_long")) ) {
				filteredList.add( alias );
			}
		}
		return filteredList;
	}

	public static Collection getSortedAliases( Collection<Alias> aliases ) {

		final AliasComparator comparator = new AliasByTypeComparator();
		Collections.sort((List<Alias>) aliases, comparator);

		return aliases;
	}

	public static boolean hasOrganism(Interactor interactor){
		if (interactor == null){
			return false;
		}

		return interactor.getOrganism() != null;
	}

	public static CrossReference getUniqueOrganismXref(Collection<CrossReference> identifiers) {
		CrossReference currentCrossReference = null;
		for (CrossReference identifier : identifiers){

			if (currentCrossReference == null){
				currentCrossReference = identifier;
			}
			if (identifier.getText() != null && currentCrossReference.getText().length() < identifier.getText().length()){
				currentCrossReference = identifier;
			}
		}

		return currentCrossReference;
	}

    public static boolean getSelectedFromMap( Map columnMap, String columnName ) {

        if ( columnMap.containsKey( columnName ) ) {
            return ( Boolean ) columnMap.get( columnName );
        }
         return false;
     }

    public static String encodeURL( String toEncode ) throws UnsupportedEncodingException {
        String s = "";
        if ( toEncode != null ) {
            s = URLEncoder.encode( toEncode, "UTF-8" );
        }
        return s;
    }

    public static String decodeURL( String toDecode ) throws UnsupportedEncodingException {
        String s = "";
        if ( toDecode != null ) {
            s = URLDecoder.decode( toDecode, "UTF-8" );
        }
        return s;
    }

}
