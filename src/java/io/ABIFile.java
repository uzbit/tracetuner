/*  1.4
 *
 *  ABIFile.java - Read an ABI sample file
 */

package com.paracel.tt.io;

import java.io.*;
import java.lang.*;

public class ABIFile extends File
{
    protected byte bytes[];
    protected int dirloc;		// offset into bytes of directory
    protected int tagcount;		// # tags in file
    protected int fileType;
    protected float scf_version_number;
    static final int ABI = 1;
    static final int SCF = 2;

    public ABIFile(String chromatFileName)
	throws FileNotFoundException, IOException
    {
        super(chromatFileName);

        FileInputStream i = new FileInputStream(this);
        int filesize = i.available();
        bytes = new byte[filesize];
        int n = i.read(bytes);
        i.close();

        if (n != filesize)
        {
            String msg = getName() + ": could read only " + n
                         + " bytes (of " + i.available() + ")";
            throw new IOException(msg);
        }

        if (bytes[0] == 'A' && bytes[1] == 'B' && bytes[2] == 'I' && bytes[3] == 'F')
            fileType = ABI;
        else if (bytes[0] == '.' && bytes[1] == 's' && bytes[2] == 'c' && bytes[3] == 'f')
            fileType = SCF;
        else
        {
            String msg = getName() + ": not an ABI or SCF file";
            throw new IOException(msg);
        }

        if (fileType == ABI)
        {
            dirloc   = get_offset(26);
            tagcount = get_offset(18);
        }
        else
        {
            String version = new String(bytes, 36, 4);
            scf_version_number = Float.valueOf(version).floatValue();
        }
    }

    public int getType()
    {
        return fileType;
    }

    public int getSize()
    {
        return bytes.length;
    }

    public int getTagCount()
    {
        return tagcount;
    }

    public int getNumCalledBases()
    {
        if (fileType == ABI)
        {
            int location = find_dir_entry("PBAS", 2);

            if (location > 0) return get_offset(location + 12);
            else return 0;
        }
        else
            return get_offset(12);
    }

    public int getNumEditedBases()
    {
        int location = find_dir_entry("PBAS", 1);

        if (location > 0) return get_offset(location + 12);
        else return 0;
    }

    public String getCalledBases()
    {
        int i;

        if (fileType == ABI)
        {
            int location = find_dir_entry("PBAS", 2);

            if (location > 0)
            {
                int base_count = get_offset(location + 12);
                int base_start = get_offset(location + 20);
                String bases = new String(bytes, base_start, base_count);

                return bases;
            }
            else
                return "";
        }
        else   /* Presumably, it's a .scf file. */
        {
            int base_count = get_offset(12);
            int base_start = get_offset(24);
            char sequence[] = new char[base_count];

            if (scf_version_number < 2.9)
                for (i = 0; i < base_count; i++)
                    sequence[i] = (char) bytes[base_start + i * 12 + 8];
            else
            {
                base_start += (base_count * 8);

                for (i = 0; i < base_count; i++)
                    sequence[i] = (char) bytes[base_start + i];
            }

            String bases = new String(sequence);

            return bases;
        }
    }

    public String getEditedBases()
    {
        int location = find_dir_entry("PBAS", 1);

        if (location > 0)
        {
            int base_count = get_offset(location + 12);
            int base_start = get_offset(location + 20);
            String bases = new String(bytes, base_start, base_count);

            return bases;
        }

        return "";
    }

    public short[] getEditedPeakLocations()
    {
        short called_locs[] = new short[0];

        if (fileType == ABI) {
            int location = find_dir_entry("PLOC", 1);

            if (location > 0) {
                int num_peaks = get_offset(location + 12);
                int peak_start = get_offset(location + 20);
                called_locs = new short[num_peaks];

                for (int i = 0; i < num_peaks; i++)
                    called_locs[i] = (short) (byteToUnsigned(bytes[peak_start + i * 2]) * 256
                                            + byteToUnsigned(bytes[peak_start + i * 2 + 1]));
                return called_locs;
            } else {
                return null;
            }
        } else {
            // SCF format file does not contain edited bases or edited base
            // peak locations data.
            return null;
        }
    }

    public short[] getCalledPeakLocations()
    {
        int i;
        short called_locs[] = new short[0];

        if (fileType == ABI)
        {
            int location = find_dir_entry("PLOC", 2);

            if (location > 0)
            {
                int num_peaks = get_offset(location + 12);
                int peak_start = get_offset(location + 20);
                called_locs = new short[num_peaks];

                for (i = 0; i < num_peaks; i++)
                    called_locs[i] = (short) (byteToUnsigned(bytes[peak_start + i * 2]) * 256
                                            + byteToUnsigned(bytes[peak_start + i * 2 + 1]));
            }
        }
        else
        {
            int num_peaks = get_offset(12);
            int peak_start = get_offset(24);
            called_locs = new short[num_peaks];

            if (scf_version_number < 2.9)
                for (i = 0; i < num_peaks; i++)
                    called_locs[i] = (short) get_offset(peak_start + i * 12);
            else
                for (i = 0; i < num_peaks; i++)
                    called_locs[i] = (short) get_offset(peak_start + i * 4);
        }

        return called_locs;
    }

    public String getBaseFromDyeIndex(int index)
    {
        int i;

        if (fileType == ABI)
        {
            int location = find_dir_entry("FWO_", 1);

            if (location > 0 && index >= 1 && index <= 4)
                return new String(bytes, location + 19 + index, 1);
            else
                return "";
        }
        else
            switch (index)   /* .scf files always have this order. */
            {
                case 1:
                    return new String("A");
                case 2:
                    return new String("C");
                case 3:
                    return new String("G");
                case 4:
                    return new String("T");
                default:
                    return "";
            }
    }

    public int[] getAnalyzedData(int dye)
    {
        int i;
        int analyzed_array[] = new int[0];

        if (fileType == ABI)
        {
            int location = find_dir_entry("DATA", dye + 8);

            if (location > 0)
            {
                int num_data = get_offset(location + 12);
                int data_start = get_offset(location + 20);
                analyzed_array = new int[num_data];

                for (i = 0; i < num_data; i++)
                    analyzed_array[i] = byteToUnsigned(bytes[data_start + i * 2]) * 256
                                      + byteToUnsigned(bytes[data_start + i * 2 + 1]);
            }
        }
        else
        {
            int num_data = get_offset(4);
            int data_start = get_offset(8);
            analyzed_array = new int[num_data];

            if (scf_version_number < 2.9)
            {
                if (get_offset(40) == 1)
                    for (i = 0; i < num_data; i++)
                        analyzed_array[i] = byteToUnsigned(bytes[data_start + i * 4 + dye - 1]);
                else
                    for (i = 0; i < num_data; i++)
                        analyzed_array[i] = byteToUnsigned(bytes[data_start + i * 8 + (dye - 1) * 2]) * 256
                                          + byteToUnsigned(bytes[data_start + 1 + i * 8 + (dye - 1) * 2]);
            }
            else
            {
                if (get_offset(40) == 1)
                {
                    for (i = 0; i < num_data; i++)
                        analyzed_array[i] = byteToUnsigned(bytes[data_start + (dye - 1) * num_data + i]);

                    delta_samples1(analyzed_array, num_data);
                }
                else
                {
					for (i = 0; i < num_data; i++)
                    {
                        analyzed_array[i] = byteToUnsigned(bytes[data_start + ((dye - 1) * num_data * 2) + (i * 2)]);
                        analyzed_array[i] *= 256;
                        analyzed_array[i] += byteToUnsigned(bytes[data_start + ((dye - 1) * num_data * 2) + (i * 2) + 1]);
                    }

                    delta_samples2(analyzed_array, num_data);
                }
            }
        }

        return analyzed_array;
    }

    public void delta_samples1(int samples[], int num_samples)
    {
        int i, p_sample = 0;

        for (i = 0; i < num_samples; i++)
        {
            samples[i] += p_sample;
            samples[i] &= 0xff;
            p_sample = samples[i];
        }

        p_sample = 0;

        for (i = 0; i < num_samples; i++)
        {
            samples[i] += p_sample;
            samples[i] &= 0xff;
            p_sample = samples[i];
        }
    }

    public void delta_samples2(int samples[], int num_samples)
    {
        int i, p_sample = 0;

        for (i = 0; i < num_samples; i++)
        {
            samples[i] += p_sample;
            samples[i] &= 0xffff;
            p_sample = samples[i];
        }

        p_sample = 0;

        for (i = 0; i < num_samples; i++)
        {
            samples[i] += p_sample;
            samples[i] &= 0xffff;
            p_sample = samples[i];
        }
    }

    public int[] getRawData(int dye)
    {
        int raw_array[] = new int[0];

        int location = find_dir_entry("DATA", dye);

        if (location > 0)
        {
            int num_data = get_offset(location + 12);
            int data_start = get_offset(location + 20);
            raw_array = new int[num_data];

            for (int i = 0; i < num_data; i++)
                raw_array[i] = byteToUnsigned(bytes[data_start + i * 2]) * 256
                             + byteToUnsigned(bytes[data_start + i * 2 + 1]);
        }

        return raw_array;
    }

    public int[] getFifthDyeData()
    {
        int fifth_dye_array[] = new int[0];

        int location = find_dir_entry("DATA", 105);

        if (location > 0)
        {
            int num_data = get_offset(location + 12);
            int data_start = get_offset(location + 20);
            fifth_dye_array = new int[num_data];

            for (int i = 0; i < num_data; i++)
                fifth_dye_array[i] = byteToUnsigned(bytes[data_start + i * 2]) * 256
                                   + byteToUnsigned(bytes[data_start + i * 2 + 1]);
        }

        return fifth_dye_array;
    }

    public String getTag(int index)
    {
        if (index < 0 || index >= tagcount) return null;

        int tagloc = dirloc + 28 * index;
        char tag[] = new char[4];
        tag[0] = (char) bytes[tagloc++];
        tag[1] = (char) bytes[tagloc++];
        tag[2] = (char) bytes[tagloc++];
        tag[3] = (char) bytes[tagloc++];

        return new String(tag);
    }

    protected int byteToUnsigned(byte b)
    {
        if (b < 0) return 256 + b;

        return b;
    }

    protected int get_offset(int index)
    {
        int offset = byteToUnsigned(bytes[index]) * 256 * 256 * 256
                   + byteToUnsigned(bytes[index + 1]) * 256 * 256
                   + byteToUnsigned(bytes[index + 2]) * 256
                   + byteToUnsigned(bytes[index + 3]);

        return offset;
    }

    protected int find_dir_entry(String tag, int id)
    {
        byte t[];

        if (tag.length() != 4) return 0;

        try { t = tag.getBytes(); }
        catch (Exception e) { return 0; }

        for (int i = 0, tagloc = dirloc; i < tagcount; i++, tagloc += 28)
            if (bytes[tagloc] == t[0] && bytes[tagloc + 1] == t[1]
            &&  bytes[tagloc + 2] == t[2] && bytes[tagloc + 3] == t[3]
            &&  bytes[tagloc + 7] == id)
            {
                return tagloc;
            }

        return 0;
    }

    public void releaseResources() {
        bytes = null;
    }

    /*
     *  The following main procedure is for debugging only.  Normally
     *  this class is called from other classes.
     */
    public static void main(String args[])
    {
        ABIFile f;
        int i, n, analyzed_data[], raw_data[];

        if (args.length < 1)
        {
            System.err.println("usage: java ABIFile <pathname>");
            System.exit(2);
        }

        try { f = new ABIFile(args[0]); }
        catch (Exception e)
        {
            System.err.println("Error opening file: " + e.getMessage());

            return;
        }

        if (f.fileType == SCF)
        {
            System.out.println(f.getName() + ": " + f.getSize() + " bytes");
            System.out.println(f.getNumCalledBases() + " called bases:");
            System.out.println(f.getCalledBases());

            analyzed_data = f.getAnalyzedData(1);

            if (analyzed_data.length > 20)
            {
                System.out.println("Analyzed data (" + analyzed_data.length + " points):");

                for (i = 0; i < 20; i++)
                    System.out.print(analyzed_data[i] + ", ");

                System.out.println();
            }

            System.exit(0);
        }

        System.out.println(f.getName() + ": " + f.getSize() + " bytes, " +
                           f.getTagCount() + " tags");
        System.out.println(f.getNumCalledBases() + " called bases:");
        System.out.println(f.getCalledBases());
        System.out.println(f.getNumEditedBases() + " edited bases:");
        System.out.println(f.getEditedBases());

        raw_data      = f.getRawData(1);
        analyzed_data = f.getAnalyzedData(1);

        if (raw_data.length > 10)
        {
            System.out.println("Raw data (" + raw_data.length + " points):");

            for (i = 0; i < 10; i++)
                System.out.print(raw_data[i] + ", ");

            System.out.println();
        }

        if (analyzed_data.length > 10)
        {
            System.out.println("Analyzed data (" + analyzed_data.length + " points):");

            for (i = 0; i < 10; i++)
                System.out.print(analyzed_data[i] + ", ");

            System.out.println();
        }
    }
}
