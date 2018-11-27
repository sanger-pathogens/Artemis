package uk.ac.sanger.artemis.components;

import static org.junit.Assert.*;
import static org.mockito.Mockito.*;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.UnknownHostException;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;
import org.mockito.Spy;

import uk.ac.sanger.artemis.Options;

/**
 * Unit test for the RunInterProSearch class.
 * 
 * @author kp11
 *
 */
public class RunInterProSearchTest
{
	private static final String TEST_WEBSITE_URL = "https:/interpro-site";
	private static final String TEST_WEBSITE_URL_REDIRECT = "https:/interpro-site/jobpage";
	private static final String TEST_SEQUENCE = "GALM";
	private static final String TEST_EXCEPTION_MESSAGE = "Oops something went wrong";
	
	@Spy
	private Options options;
	
	@Mock 
	private HttpURLConnection conn;
	@Mock 
	private OutputStream outStream;
	@Mock
	private InputStream  inStream;
	
	@Spy
	private RunInterProSearch ipSearch;
	
	@Before
	public void setUp() throws Exception {

		MockitoAnnotations.initMocks(this);
	}
	
	@Test
	public void testGetUrlFromOptions() 
	{
		// Just check the InterPro search URL is defined in the options file.
		// Don't bother to check the actual String.
		
		assertNotNull(options.getProperty("interpro_search_url"));
		assertTrue(options.getProperty("interpro_search_url").length() > 0);
	}
	
	@Test
	public void testRunInterProSearchConstructors() 
	{
		RunInterProSearch interproSearch1 = new RunInterProSearch(TEST_SEQUENCE);
		assertEquals(TEST_SEQUENCE, interproSearch1.getSequence());
		assertTrue(interproSearch1.isDaemon());
		
		RunInterProSearch interproSearch2 = new RunInterProSearch(TEST_SEQUENCE, TEST_WEBSITE_URL);
		assertEquals(TEST_SEQUENCE, interproSearch2.getSequence());
		assertEquals(TEST_WEBSITE_URL, interproSearch2.getSearchURL());
	}
	
	@Test
	public void testStart1() throws Exception 
	{
		// Given
		ipSearch.setSequence(TEST_SEQUENCE);
		
		// When
		when( options.getProperty("interpro_search_url") ).thenReturn(TEST_WEBSITE_URL);
		when( ipSearch.getConnection() ).thenReturn(conn);
		when( conn.getOutputStream() ).thenReturn(outStream);
		when( conn.getInputStream() ).thenReturn(inStream);
		when( conn.getHeaderField(anyString()) ).thenReturn(TEST_WEBSITE_URL_REDIRECT);
		when( conn.getResponseCode() ).thenReturn(HttpURLConnection.HTTP_OK);
		doNothing().when(ipSearch).displayURL(any()); // So we don't pop up a browser window
		
		ipSearch.run();
		
		// Then
		verify(ipSearch, times(1)).getConnection();
		verify(conn, times(1)).getOutputStream();
		verify(conn, times(1)).getInputStream();
		verify(conn, times(1)).setDoOutput(eq(true));
		verify(conn, times(1)).setInstanceFollowRedirects(eq(false));
		verify(conn, times(1)).getResponseCode();
		verify(ipSearch, times(1)).displayURL(TEST_WEBSITE_URL_REDIRECT);
		
		// Verify clean up.
		verify(outStream, times(1)).close();
		verify(inStream, times(1)).close();
	}
	
	@Test
	public void testStart2() throws Exception 
	{
		// Given
		ipSearch.setSequence(TEST_SEQUENCE);
		
		// When
		when( options.getProperty("interpro_search_url") ).thenReturn(TEST_WEBSITE_URL);
		when( ipSearch.getConnection() ).thenReturn(conn);
		when( conn.getOutputStream() ).thenReturn(outStream);
		when( conn.getInputStream() ).thenReturn(inStream);
		when( conn.getHeaderField(anyString()) ).thenReturn(TEST_WEBSITE_URL_REDIRECT);
		when( conn.getResponseCode() ).thenReturn(HttpURLConnection.HTTP_MOVED_TEMP); // Different response code
		doNothing().when(ipSearch).displayURL(any()); // So we don't pop up a browser window
		
		ipSearch.run();
		
		// Then
		verify(ipSearch, times(1)).getConnection();
		verify(conn, times(1)).getOutputStream();
		verify(conn, times(1)).getInputStream();
		verify(conn, times(1)).setDoOutput(eq(true));
		verify(conn, times(1)).setInstanceFollowRedirects(eq(false));
		verify(conn, times(1)).getResponseCode();
		verify(ipSearch, times(1)).displayURL(TEST_WEBSITE_URL_REDIRECT);
		
		// Verify clean up.
		verify(outStream, times(1)).close();
		verify(inStream, times(1)).close();		
	}
	
	@Test
	public void testStartWithUnknownHostException() throws Exception 
	{
		// Given
		ipSearch.setSequence(TEST_SEQUENCE);
		
		// When
		when( options.getProperty("interpro_search_url") ).thenReturn(TEST_WEBSITE_URL);
		when( ipSearch.getConnection() ).thenThrow(new UnknownHostException(TEST_EXCEPTION_MESSAGE));
		when( conn.getOutputStream() ).thenReturn(outStream);
		when( conn.getInputStream() ).thenReturn(inStream);
		when( conn.getHeaderField(anyString()) ).thenReturn(TEST_WEBSITE_URL_REDIRECT);
		when( conn.getResponseCode() ).thenReturn(HttpURLConnection.HTTP_OK);
		doNothing().when(ipSearch).displayURL(any()); 
		doNothing().when(ipSearch).showError(any()); // mock out error display
		
		ipSearch.run();
		
		// Then
		verify(ipSearch, times(1)).showError(anyString());
		
		// Verify clean up.
		verify(outStream, times(0)).close();
		verify(inStream, times(0)).close();
		
	}
	
	@Test
	public void testStartWithGeneralException() throws Exception 
	{
		// Given
		ipSearch.setSequence(TEST_SEQUENCE);
		
		// When
		when( options.getProperty("interpro_search_url") ).thenReturn(TEST_WEBSITE_URL);
		when( ipSearch.getConnection() ).thenThrow(new IOException(TEST_EXCEPTION_MESSAGE));
		when( conn.getOutputStream() ).thenReturn(outStream);
		when( conn.getInputStream() ).thenReturn(inStream);
		when( conn.getHeaderField(anyString()) ).thenReturn(TEST_WEBSITE_URL_REDIRECT);
		when( conn.getResponseCode() ).thenReturn(HttpURLConnection.HTTP_OK);
		doNothing().when(ipSearch).displayURL(any());
		doNothing().when(ipSearch).showError(any()); // mock out error display
		
		ipSearch.run();
		
		// Then
		verify(ipSearch, times(1)).showError(anyString());
		
		// Verify clean up.
		verify(outStream, times(0)).close();
		verify(inStream, times(0)).close();
		
	}
	
	@Test
	public void testStartWithBadResponseCode() throws Exception 
	{
		// Given
		ipSearch.setSequence(TEST_SEQUENCE);
		
		// When
		when( options.getProperty("interpro_search_url") ).thenReturn(TEST_WEBSITE_URL);
		when( ipSearch.getConnection() ).thenReturn(conn);
		when( conn.getOutputStream() ).thenReturn(outStream);
		when( conn.getInputStream() ).thenReturn(inStream);
		when( conn.getHeaderField(anyString()) ).thenReturn(TEST_WEBSITE_URL_REDIRECT);
		when( conn.getResponseCode() ).thenReturn(HttpURLConnection.HTTP_FORBIDDEN);
		doNothing().when(ipSearch).displayURL(any());
		doNothing().when(ipSearch).showError(any()); // mock out error display
		
		ipSearch.run();
		
		// Then
		verify(ipSearch, times(1)).showError(anyString());
		
		// Verify clean up.
		verify(outStream, times(1)).close();
		verify(inStream, times(1)).close();
		
	}
}
