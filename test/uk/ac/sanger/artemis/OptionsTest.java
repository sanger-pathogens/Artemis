/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package test.uk.ac.sanger.artemis;

import org.junit.Test;
import static org.junit.Assert.*;
import uk.ac.sanger.artemis.*;

/**
 *
 * @author luj
 */
public class OptionsTest {
    
    public OptionsTest() {
    }

    /**
     * Test of getPropertyTruthValue method, of class Options.
     */
    @Test
    public void testcanViewCustomAnnotation() {
        System.setProperty("view_Custom_Annotation", "true");
        Options op = new Options();
        assertEquals(true, op.canViewCustomAnnotation());
    }

    
}