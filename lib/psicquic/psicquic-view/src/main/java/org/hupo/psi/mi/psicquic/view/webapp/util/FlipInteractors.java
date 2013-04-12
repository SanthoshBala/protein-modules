package org.hupo.psi.mi.psicquic.view.webapp.util;

import psidev.psi.mi.tab.model.Alias;
import psidev.psi.mi.tab.model.BinaryInteraction;
import psidev.psi.mi.tab.model.CrossReference;
import psidev.psi.mi.tab.model.Interactor;

import java.util.Collection;

/**
 * Created with IntelliJ IDEA.
 * User: ntoro
 * Date: 09/10/2012
 * Time: 11:37
 * To change this template use File | Settings | File Templates.
 */
public class FlipInteractors {

	public void flipIfNecessary(BinaryInteraction interaction, String searchQuery) {

		final Interactor interactorA = interaction.getInteractorA();
		final Interactor interactorB = interaction.getInteractorB();

		if (interactorA != null && interactorB != null) {

			final boolean matchesA = matchesQuery(interactorA, searchQuery);
			final boolean matchesB = matchesQuery(interactorB, searchQuery);

			if (matchesA && !matchesB) {
				// nothing
			} else if (!matchesA && matchesB) {
				interaction.flip();
			} else {
				String id1 = (!interactorA.getIdentifiers().isEmpty() ? interactorA.getIdentifiers().iterator().next().getIdentifier() : "");
				String id2 = (!interactorB.getIdentifiers().isEmpty() ? interactorB.getIdentifiers().iterator().next().getIdentifier() : "");

				int comp = id1.compareTo(id2);
				if (comp > 0) {
					interaction.flip();
				}
			}
		}
	}

	private boolean matchesQuery(Interactor interactor, String searchQuery) {

		if (interactor != null) {
			String queries[] = searchQuery.split(" ");

			for (String query : queries) {
				if ("NOT".equalsIgnoreCase(query) ||
						"AND".equalsIgnoreCase(query) ||
						"OR".equalsIgnoreCase(query)) {
					continue;
				}

				if (matchesQueryAliases(query, interactor.getAliases())) {
					return true;
				} else if (matchesQueryXrefs(query, interactor.getIdentifiers())) {
					return true;
				} else if (matchesQueryXrefs(query, interactor.getAlternativeIdentifiers())) {
					return true;
				} else if (matchesQueryXrefs(query, interactor.getXrefs())) {
					return true;
				}
			}
		}
		return false;
	}

	private boolean matchesQueryAliases(String query, Collection<Alias> aliases) {
		for (Alias alias : aliases) {
			if (query.toLowerCase().toLowerCase().endsWith(alias.getName().toLowerCase())) {
				return true;
			}
		}
		return false;
	}

	private boolean matchesQueryXrefs(String query, Collection<CrossReference> xrefs) {
		for (CrossReference xref : xrefs) {
			if (query.toLowerCase().endsWith(xref.getIdentifier().toLowerCase())) {
				return true;
			}
		}
		return false;
	}
}
