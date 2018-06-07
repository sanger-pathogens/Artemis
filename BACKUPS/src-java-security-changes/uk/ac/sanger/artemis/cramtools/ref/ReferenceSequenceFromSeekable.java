/*******************************************************************************
 * Copyright 2013 EMBL-EBI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

package uk.ac.sanger.artemis.cramtools.ref;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.nio.ByteBuffer;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.MatchResult;


/**
 * This is to compensate the lack of ReferenceSequenceFile working with a URL.
 * 
 * @author vadim
 * 
 */
class ReferenceSequenceFromSeekable {

	private static final int BUFFER_SIZE = 1024 * 100;
	private SeekableStream s;
	private Map<String, FastaSequenceIndexEntry> index;

	private ReferenceSequenceFromSeekable(SeekableStream s, Map<String, FastaSequenceIndexEntry> index) {
		this.s = s;
		this.index = index;
	}

	public byte[] getSubsequenceAt(String contig, long start, long stop) {
		if (start > stop + 1)
			throw new SAMException(
					String.format("Malformed query; start point %d lies after end point %d", start, stop));

		FastaSequenceIndexEntry indexEntry = index.get(contig);

		if (stop > indexEntry.getSize())
			throw new SAMException("Query asks for data past end of contig");

		int length = (int) (stop - start + 1);

		byte[] target = new byte[length];
		ByteBuffer targetBuffer = ByteBuffer.wrap(target);

		final int basesPerLine = indexEntry.getBasesPerLine();
		final int bytesPerLine = indexEntry.getBytesPerLine();
		final int terminatorLength = bytesPerLine - basesPerLine;

		long startOffset = ((start - 1) / basesPerLine) * bytesPerLine + (start - 1) % basesPerLine;

		// Allocate a 128K buffer for reading in sequence data.
		ByteBuffer channelBuffer = ByteBuffer.allocate(BUFFER_SIZE);

		while (targetBuffer.position() < length) {
			// If the bufferOffset is currently within the eol characters in the
			// string, push the bufferOffset forward to the next printable
			// character.
			startOffset += Math.max((int) (startOffset % bytesPerLine - basesPerLine + 1), 0);

			try {
				s.seek(indexEntry.getLocation() + startOffset);
				startOffset += Utils.readInto(channelBuffer, s);
			} catch (IOException ex) {
				throw new SAMException("Unable to load " + contig + "(" + start + ", " + stop + ") from "
						+ s.getSource());
			}

			// Reset the buffer for outbound transfers.
			channelBuffer.flip();

			// Calculate the size of the next run of bases based on the contents
			// we've already retrieved.
			final int positionInContig = (int) start - 1 + targetBuffer.position();
			final int nextBaseSpan = Math.min(basesPerLine - positionInContig % basesPerLine,
					length - targetBuffer.position());
			// Cap the bytes to transfer by limiting the nextBaseSpan to the
			// size of the channel buffer.
			int bytesToTransfer = Math.min(nextBaseSpan, channelBuffer.capacity());

			channelBuffer.limit(channelBuffer.position() + bytesToTransfer);

			while (channelBuffer.hasRemaining()) {
				targetBuffer.put(channelBuffer);

				bytesToTransfer = Math.min(basesPerLine, length - targetBuffer.position());
				channelBuffer.limit(Math.min(channelBuffer.position() + bytesToTransfer + terminatorLength,
						channelBuffer.capacity()));
				channelBuffer.position(Math.min(channelBuffer.position() + terminatorLength, channelBuffer.capacity()));
			}

			// Reset the buffer for inbound transfers.
			channelBuffer.flip();
		}

		return target;
	}

	public static ReferenceSequenceFromSeekable fromString(String urlOrPathString) {
		try {
			SeekableStream fastaSS = SeekableStreamFactory.getInstance().getStreamFor(urlOrPathString);
			SeekableStream indexSS = SeekableStreamFactory.getInstance().getStreamFor(urlOrPathString + ".fai");
			return new ReferenceSequenceFromSeekable(fastaSS, buildIndex(indexSS));
		} catch (MalformedURLException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private static Map<String, FastaSequenceIndexEntry> buildIndex(InputStream is) {
		Scanner scanner = new Scanner(is);

		int sequenceIndex = 0;
		Map<String, FastaSequenceIndexEntry> index = new HashMap<String, FastaSequenceIndexEntry>();
		while (scanner.hasNext()) {
			// Tokenize and validate the index line.
			String result = scanner.findInLine("(.+)\\t+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)");
			if (result == null)
				throw new RuntimeException("Found invalid line in index file:" + scanner.nextLine());
			MatchResult tokens = scanner.match();
			if (tokens.groupCount() != 5)
				throw new RuntimeException("Found invalid line in index file:" + scanner.nextLine());

			// Skip past the line separator
			scanner.nextLine();

			// Parse the index line.
			String contig = tokens.group(1);
			long size = Long.valueOf(tokens.group(2));
			long location = Long.valueOf(tokens.group(3));
			int basesPerLine = Integer.valueOf(tokens.group(4));
			int bytesPerLine = Integer.valueOf(tokens.group(5));

			contig = SAMSequenceRecord.truncateSequenceName(contig);
			// Build sequence structure
			index.put(contig, new FastaSequenceIndexEntry(contig, location, size, basesPerLine, bytesPerLine,
					sequenceIndex++));
		}
		scanner.close();
		return index;
	}

	private static class FastaSequenceIndexEntry {
		private String contig;
		private long location;
		private long size;
		private int basesPerLine;
		private int bytesPerLine;
		private final int sequenceIndex;

		public FastaSequenceIndexEntry(String contig, long location, long size, int basesPerLine, int bytesPerLine,
				int sequenceIndex) {
			this.contig = contig;
			this.location = location;
			this.size = size;
			this.basesPerLine = basesPerLine;
			this.bytesPerLine = bytesPerLine;
			this.sequenceIndex = sequenceIndex;
		}

		public String getContig() {
			return contig;
		}

		public long getLocation() {
			return location;
		}

		public long getSize() {
			return size;
		}

		public int getBasesPerLine() {
			return basesPerLine;
		}

		public int getBytesPerLine() {
			return bytesPerLine;
		}

		public int getSequenceIndex() {
			return sequenceIndex;
		}
	}
}
