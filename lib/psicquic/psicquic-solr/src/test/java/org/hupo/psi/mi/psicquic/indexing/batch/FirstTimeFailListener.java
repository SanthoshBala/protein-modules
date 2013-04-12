/**
 * Copyright 2009 The European Bioinformatics Institute, and others.
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
package org.hupo.psi.mi.psicquic.indexing.batch;

import org.springframework.batch.core.ItemReadListener;

/**
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 */
public class FirstTimeFailListener implements ItemReadListener {

    private boolean secondTime = false;
    private int count = 0;

    public void beforeRead() {
        if (count == 2 && !secondTime) {
            secondTime = true;
            throw new RuntimeException("Breaking thing");
        }
    }

    public void afterRead(Object item) {
        count++;
    }

    public void onReadError(Exception ex) {

    }
}
