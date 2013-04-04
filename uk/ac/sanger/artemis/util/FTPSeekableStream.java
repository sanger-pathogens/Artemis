package uk.ac.sanger.artemis.util;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.SocketException;
import java.net.SocketTimeoutException;
import java.net.URL;

import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPConnectionClosedException;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.commons.net.ftp.FTPReply;
import org.apache.log4j.Logger;

import net.sf.samtools.seekablestream.SeekableStream;

/**
 * Written independently to, but bugfixed by looking at the Savant
 * SeekableFTPStream.
 * 
 * @author gv1
 * 
 */
public class FTPSeekableStream extends SeekableStream {

    private static final Logger logger = Logger
            .getLogger(FTPSeekableStream.class);

    private static final String defaultUser = "anonymous";
    private static final String defaultPassword = "";

    private URL url;
    private String host;

    private String user;
    private String password;

    private String remoteFilePath;
    private String remoteFileName;

    private FTPClient _client;

    private long position = 0;
    private long length = -1;

    private File tmpFolder;
    private File index;

    public FTPSeekableStream(URL url) throws SocketException, IOException {
        this(url, defaultUser, defaultPassword);
    }

    public FTPSeekableStream(URL url, String user, String password)
            throws SocketException, IOException {

        this.url = url;
        this.user = user;
        this.password = password;

        host = url.getHost();
        remoteFilePath = url.getPath();

        String[] split = remoteFilePath.split("/");
        if (split.length > 0) {
            remoteFileName = split[split.length - 1];
        } else {
            remoteFileName = remoteFilePath;
        }

        logger.info(String.format("Setup a stream for %s %s %s %s", host,
                remoteFilePath, this.user, this.password));

    }

    private FTPClient getClient() throws SocketException, IOException {

        if (_client == null) {

            FTPClient client = new FTPClient();

            client.connect(host);
            client.login(this.user, this.password);

            logger.debug(client.getReplyString());

            client.setFileType(FTP.BINARY_FILE_TYPE);
            client.enterLocalPassiveMode();
            client.setSoTimeout(10000);

            int reply = client.getReplyCode();
            logger.info(reply);

            if (!FTPReply.isPositiveCompletion(reply)) {
                close();
                throw new IOException("FTP server refused connection.");
            }

            _client = client;
        }
        return _client;
    }

    @Override
    public long length() {
        logger.info("length " + length);

        if (length != -1) {
            return length;
        }

        try {
            for (FTPFile f : getClient().listFiles(remoteFilePath)) {
                if (f.getName().equals(remoteFilePath)) {
                    length = f.getSize();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return length;

    }

    @Override
    public int read(byte[] bytes, int offset, int length)
            throws IOException {

        InputStream in = initStream();

        if (in == null) {
            throw new IOException("Could not get stream");
        }

        int i = 0;

        while (i < length) {
            int bytesRead = in.read(bytes, offset + i, length - i);
            if (bytesRead < 0) {
                if (i == 0) {
                    return -1;
                } else {
                    break;
                }
            }

            i += bytesRead;
        }

        finishStream(in);

        position += i;

        return i;
    }

    @Override
    public void close() throws IOException {

        if (_client != null) {

            try {
                _client.completePendingCommand();
            } catch (IOException e) {
                logger.error(e);
            }
            try {
                _client.logout();
            } catch (IOException e) {
                logger.error(e);
            }

            _client.disconnect();
            _client = null;
        }
    }

    @Override
    public boolean eof() throws IOException {
        if (position >= length()) {
            return true;
        }
        return false;
    }

    @Override
    public String getSource() {
        return url.toString();
    }

    @Override
    public void seek(long position) throws IOException {
        logger.info("seek " + position);
        this.position = position;

    }

    @Override
    public int read() throws IOException {
        logger.info("read");

        InputStream in = initStream();

        int read = in.read();
        position++;

        finishStream(in);

        return read;

    }

    private InputStream initStream() throws SocketException, IOException {
        FTPClient client = getClient();
        client.setRestartOffset(position);
        InputStream in = client.retrieveFileStream(remoteFilePath);
        return in;
    }

    private void finishStream(InputStream in) throws IOException {
        in.close();
        try {
            getClient().completePendingCommand();
        } catch (FTPConnectionClosedException suppressed) {
        } catch (SocketTimeoutException stx) {
            close();
        }

    }

    public File getTmpFolder() {
        if (tmpFolder == null) {
            tmpFolder = new File("/tmp/");
        }
        return tmpFolder;
    }

    public void setTmpFolder(File tmpFolder) throws IOException {
        if (!tmpFolder.isDirectory()) {
            throw new IOException("File " + tmpFolder.getName()
                    + " is not a folder");
        }
        this.tmpFolder = tmpFolder;
    }

    public File getIndexFile() throws IOException {

        if (index == null) {

            String indexFileName = remoteFileName + ".bai";
            String localPath = getTmpFolder().getAbsolutePath() + "/"
                    + indexFileName;

            File localFile = new File(localPath);

            if (!localFile.isFile()) {

                String remotePath = remoteFilePath + ".bai";
                logger.info(String.format("Downloading from %s to %s",
                        remotePath, localPath));

                FTPClient client = getClient();
                client.setRestartOffset(0);
                

                InputStream in = null;
                FileOutputStream out = null;
                
                try {

                    in = client.retrieveFileStream(remotePath);
                    out = new FileOutputStream(localPath);

                    byte[] buffer = new byte[1024];

                    int len;
                    int total = 0;

                    while ((len = in.read(buffer)) > 0) {
                        total += len;
                        out.write(buffer, 0, len);
                    }

                    logger.info("Index Size in bytes : " + total);

                    client.completePendingCommand();

                    index = new File(localPath);

                    if (!index.isFile()) {
                        throw new IOException(
                                "Could not save the index file locally");
                    }

                    logger.info("Saved " + index.getAbsolutePath());

                } finally {
                    
                    if (in != null) {
                        in.close();
                    }
                    if (out != null) {
                        out.close();
                    }
                    
                }

            } else {
                logger.info("Using cached index " + localFile.getAbsolutePath());
                index = localFile;
            }

            logger.info("File size " + index.length());

        }

        return index;

    }

    @Override
    public long position() throws IOException
    {
      // TODO Auto-generated method stub
      return 0;
    }

}
