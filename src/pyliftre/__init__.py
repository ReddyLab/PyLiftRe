import gzip

from bx.intervals.intersection import Interval, Intersecter


class Mapper:
    def __init__(self, chain_file: str):
        self.chain = self._read_chain_file(chain_file)

    def _update_chromID(self, c_temp, c_target, chr_style="a") -> str:
        """
        Update chromsome ID styles from 'c_target' to 'c_temp'.

        Parameters
        ----------
        c_temp : str
            Template of chromsome ID

        c_target : str
            Chromosome ID that need to be updated

        chr_style : str, optional
            Chromosome ID style. Must be one of ['a', 's', 'l'], where
            'a' : as-is. The chromosome ID of the target is written to the output file in the same \
                style of the template chromosome ID (which comes from the query). This is a \
                per-function call effect, and, thus, it is applied individually to each \
                query-ID/target-ID pair presented in any given input record. The output \
                file may have mixed styles if the input file has mixed styles.
            's' : short ID, such as "1", "2", "X.
            'l' : long ID, such as "chr1", "chr2", "chrX.
            'n' : no-change. The chromosome ID is left completely unchanged. \
                c_temp is completely ignored.

        Returns
        --------
        Updated chromosome ID

        Examples
        --------
        >>> update_chromID('chrX',1)
        'chr1'
        >>> update_chromID('1','chrY')
        'Y'
        """
        c_temp = str(c_temp)
        c_target = str(c_target)

        match chr_style:
            case "n":  # no change
                return c_target
            case "s":  # short style
                return c_target[3:] if c_target.startswith("chr") else c_target
            case "l":  # long style
                return c_target if c_target.startswith("chr") else f"chr{c_target}"
            case _:  # as-is
                if c_temp.startswith("chr"):
                    return c_target if c_target.startswith("chr") else f"chr{c_target}"
                else:
                    return c_target[3:] if c_target.startswith("chr") else c_target

    def _intersect_regions(self, lst1, lst2):
        """
        Return intersection of two bed regions.

        Parameters
        ----------
        lst1 : list
            The 1st genomic region. List of chrom, start, end.
            Example: ['chr1',10, 100]

        lst2 : list
            The 2nd genomic region. List of chrom, start, end.
            Example: ['chr1',50, 120]

        Examples
        --------
        >>> intersect_regions(['chr1',10, 100],['chr1',50, 120])
        ('chr1', 50, 100)
        >>> intersect_regions(['chr1',10, 100],['chr1',20, 30])
        ('chr1', 20, 30)

        """
        (chr1, st1, end1), (chr2, st2, end2) = lst1, lst2

        if int(st1) > int(end1) or int(st2) > int(end2):
            raise ValueError(f"Start cannot be larger than end : {lst1}, {lst2}")

        if chr1 != chr2:
            return None

        if int(st1) > int(end2) or int(st2) > int(end1):
            return None

        return (chr1, max(st1, st2), min(end1, end2))

    def _read_chain_file(self, chain_file: str):
        """
        Read chain file.

        Parameters
        ----------
        chain_file : file
            Chain format file. Input chain_file could be either plain text, \
            compressed file (".gz",".Z", ".z", ".bz", ".bz2", ".bzip2"), or a \
            URL pointing to the chain file ("http://","https://", "ftp://"). \
            If url was used, chain file must be plain text.

        print_table : bool, optional
            Print mappings in human readable table.

        Returns
        -------
        maps : dict
            Dictionary with source chrom name as key, IntervalTree object as \
            value. An IntervalTree contains many intervals. An interval is a \
            start and end position and a value. eg. Interval(11, 12, strand="-",\
            value = "abc")
        """

        maps = {}

        c_open = gzip.open if chain_file.endswith(".gz") else open
        with c_open(chain_file, "rt") as chains:
            for line in chains:
                line = line.strip()
                if not line or line.startswith(("#", " ")):
                    continue

                fields = line.split()

                if fields[0] == "chain" and len(fields) in [12, 13]:
                    source_chrom = fields[2]  # E.g. chrY
                    source_strand = fields[4]  # Must be +
                    if source_strand != "+":
                        raise ValueError(f"Source strand in a chain file must be +. ({line})")
                    source_start = int(fields[5])  # Start of source region

                    target_chrom = fields[7]  # E.g. chr5
                    target_size = int(fields[8])  # Full length of the chromosome
                    target_strand = fields[9]  # + or -
                    target_start = int(fields[10])

                    if target_strand not in ["+", "-"]:
                        raise ValueError(f"Target strand must be - or +. ({line})")
                    # chain_id = None if len(fields) == 12 else fields[12]
                    if source_chrom not in maps:
                        maps[source_chrom] = Intersecter()

                    sfrom, tfrom = source_start, target_start

                # Now read the alignment chain from the file and store it as a list
                # (source_from, source_to) -> (target_from, target_to)
                elif len(fields) == 3:
                    size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])

                    if target_strand == "+":
                        maps[source_chrom].add_interval(
                            Interval(
                                sfrom,
                                sfrom + size,
                                (target_chrom, tfrom, tfrom + size, target_strand),
                            )
                        )
                    elif target_strand == "-":
                        maps[source_chrom].add_interval(
                            Interval(
                                sfrom,
                                sfrom + size,
                                (
                                    target_chrom,
                                    target_size - (tfrom + size),
                                    target_size - tfrom,
                                    target_strand,
                                ),
                            )
                        )

                    sfrom += size + sgap
                    tfrom += size + tgap

                elif len(fields) == 1:
                    size = int(fields[0])

                    if target_strand == "+":
                        maps[source_chrom].add_interval(
                            Interval(
                                sfrom,
                                sfrom + size,
                                (target_chrom, tfrom, tfrom + size, target_strand),
                            )
                        )
                    elif target_strand == "-":
                        maps[source_chrom].add_interval(
                            Interval(
                                sfrom,
                                sfrom + size,
                                (
                                    target_chrom,
                                    target_size - (tfrom + size),
                                    target_size - tfrom,
                                    target_strand,
                                ),
                            )
                        )
                else:
                    raise ValueError(f"Invalid chain format. ({line})")

        return maps

    def map_coordinates(self, q_chr: str, q_start: int, q_end: int, q_strand: str = "+", chrom_style: str = "a"):
        """
        Map coordinates from source (i.e. original) assembly to target \
        (i.e. new) assembly.

        Parameters
        ----------
        mapping : dict
            Dictionary with source chrom name as key, IntervalTree object as value.

        q_chr : str
            Chromosome ID of query interval

        q_start : int
            Start position of query interval.

        q_end : int
            End position of query interval.

        q_strand : str
            Strand of query interval.

        chrom_style : str, optional
            Chromosome ID style. Must be one of ['a', 's', 'l', 'n'], where
            'a' : as-is. The chromosome ID of the output is in the same style
                of the input query chromosome ID. This is applied invidivually
                to each query-ID/target-ID pair (as found in any given input
                to this function). If this function is called multiple times
                with different query-ID styles, the aggregated outputs will
                also have mixed styles.
            's' : short ID, such as "1", "2", "X.
            'l' : long ID, such as "chr1", "chr2", "chrX.
            'n' : no-change. do not modify the Chromosome ID in any way.
        """

        matches = []
        complement = {"+": "-", "-": "+"}

        if q_chr in self.chain:
            mapping = self.chain[q_chr]
        elif q_chr[3:] in self.chain:
            mapping = self.chain[q_chr[3:]]
        elif f"chr{q_chr}" in self.chain:
            mapping = self.chain[f"chr{q_chr}"]
        else:
            return None

        targets = mapping.find(q_start, q_end)

        if len(targets) == 0:
            return None
        else:
            for t in targets:
                s_start = t.start
                s_end = t.end
                t_chrom = t.value[0]
                t_chrom = self._update_chromID(q_chr, t_chrom, chr_style=chrom_style)
                t_start = t.value[1]
                t_end = t.value[2]
                t_strand = t.value[3]

                region_intersection = self._intersect_regions((q_chr, q_start, q_end), (q_chr, s_start, s_end))
                if region_intersection is None:
                    continue

                (chr, real_start, real_end) = region_intersection

                l_offset = abs(real_start - s_start)
                size = abs(real_end - real_start)
                matches.append((chr, real_start, real_end, q_strand))

                if t_strand == "+":
                    i_start = t_start + l_offset
                    if q_strand == "+":
                        matches.append((t_chrom, i_start, i_start + size, t_strand))
                    else:
                        matches.append((t_chrom, i_start, i_start + size, complement[t_strand]))
                elif t_strand == "-":
                    i_start = t_end - l_offset - size
                    if q_strand == "+":
                        matches.append((t_chrom, i_start, i_start + size, t_strand))
                    else:
                        matches.append((t_chrom, i_start, i_start + size, complement[t_strand]))
                else:
                    raise ValueError(f"Unknown strand: {q_strand}. Can only be '+' or '-'.")

        return matches
