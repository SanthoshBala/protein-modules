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
package org.hupo.psi.mi.psicquic.view.webapp.controller.search;


/**
 * Utility class to handle userQuery object
 *
 * @author Prem Anand (prem@ebi.ac.uk)
 * @version $Id$
 * @since 2.0.1-SNAPSHOT
 */
public class UserQueryUtils {

	public static void escapeIfNecessary(String query, StringBuffer queryString) {
		if (query.equals(UserQuery.STAR_QUERY)){
			queryString.append(query);
			return;
		}

		// range query, do nothing
		if (query.startsWith("[") && query.endsWith("]")){
			queryString.append(query);
		}
		else if (query.contains(" ") ||
				query.contains(":") ||
				query.contains("(") ||
				query.contains(")") ||
				query.contains("-") ||
				query.contains("+")) {

			// deal with wild search
			if (query.contains("*")){
				queryString.append(query.toLowerCase()
						.replaceAll(" ", "\\\\ ")
						.replaceAll(":", "\\\\:")
						.replaceAll("\\(", "\\\\(")
						.replaceAll("\\)", "\\\\)")
						.replaceAll("-", "\\\\-")
						.replaceAll("\\+", "\\\\+"));
			}
			else if (query.contains("(") ||
					query.contains(")") ||
					query.contains("-") ||
					query.contains("+")){
				queryString.append(query);
			}
			else {
				queryString.append("\""+query+"\"");
			}
		}
		else if (query.contains("*")){
			queryString.append(query.toLowerCase());
		}
		else {
			queryString.append(query);
		}
	}

}
