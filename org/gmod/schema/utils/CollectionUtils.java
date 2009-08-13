package org.gmod.schema.utils;

import java.util.Collection;
import java.util.HashSet;

/**
 * Static utility class for miscellaneous collection handling
 * 
 * @author art
 */
public class CollectionUtils {

    /**
     * 
     * 
     * @param <T> 
     * @param collection 
     * @return the original collection, or an empty collection
     */
    public static <T> Collection<T> safeGetter(Collection<T> collection) {
        if (collection != null) {
            return collection;
        }
        return new HashSet<T>(0);
    }
    
}
