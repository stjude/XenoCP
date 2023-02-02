#! /usr/bin/env python3

import sys
import argparse

import raptr


def main():
    """ Handles arguments and invoke driver function.
    """
    head_description = "Add a bam pair with qualifier 'xenocp' for a sample."
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    # Create a parent parser with common arguments for every subparser
    parent_parser = argparse.ArgumentParser(description="Add a bam pair with qualifier 'xenocp' for a sample.", add_help=False)
    parent_parser.add_argument('-v', '--verbose', help='Enable verbose mode', action='store_true')
    parent_parser.add_argument('-d', '--debug', help='Enable debug mode', action='store_true')
    parent_parser.add_argument('--dry-run', help='Enable dry-run mode', action='store_true')

    subparsers = parser.add_subparsers(title='Subcommands', help='Valid subcommands.', dest='subparser_name')

    # Create a subparser for arguments of sample, target, project names
    subparser_name = subparsers.add_parser('name', help='Use a sample, target, project name.', parents=[parent_parser])
    subparser_name.add_argument('-s', '--sample', metavar='', help='Sample formal_name')
    subparser_name.add_argument('-t', '--target', metavar='', help='Target name')
    subparser_name.add_argument('-p', '--project', metavar='', help='Project/subproject')

    # Create a subparser for argument of bam_id
    subparser_id = subparsers.add_parser('id', help='Use a bam_id.', parents=[parent_parser])
    subparser_id.add_argument('-i', '--id', metavar='', help='A bam_id')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if len(sys.argv) == 2:
        if args.subparser_name == 'name':
            subparser_name.print_help()
            sys.exit(0)
        elif args.subparser_name == 'id':
            subparser_id.print_help()
            sys.exit(0)

    add_bam_pair_xenocp(args)


def add_bam_pair_xenocp(args):
    """ Driver function of this script.
    Args:
        args (obj): argparse object
    Returns:
        None
    """
    rp = raptr.Raptr()
    if 'sample' in args and 'target' in args and 'project' in args:
        tokens = args.project.split('/')
        project = tokens[0]
        subproject = tokens[1]
        add_bam_pair_stp(rp, args.sample, args.target, project, subproject, dry_run=args.dry_run)
    elif 'id' in args:
        add_bam_pair_bam_id(rp, args.id, dry_run=args.dry_run)
    else:
        print('Error - miss required arguments', file=sys.stderr)

    if args.dry_run:
        rp.conn_rollback()
    else:
        rp.conn_commit()


def add_bam_pair_stp(raptr, sample, target, project, subproject, **kwargs):
    """ Add a bam pair for a given stp.
    Args:
        raptr (obj): a Raptr instance
        sample (str): a formal_name
        target (str): a target_name
        project (str): project name
        subproject (str): subproject name
        kwargs (dict): keyword argument
    Returns:
        None
    """
    # Get stp ID
    query = """select sample_target_project_id from sample_target_project_view where 
            formal_name = %s and target_name = %s and project_name = %s and subproject = %s;"""
    stp_id = raptr.fetch_item_or_fail(query, (sample, target, project, subproject))

    # Get read groups associated with old bam_tpl
    query = """select bam_tpl_id, bam_id from bam_and_tpl where sample_target_project_id = %s
            and bam_status = 'Normal' and legacy = False"""
    (bam_tpl_id, bam_id) = raptr.fetch_row_or_fail(query, (stp_id,))
    query = """select read_group_id from bam_tpl_read_group where bam_tpl_id = %s;"""
    rg_ids = raptr.fetch_col_or_fail(query, (bam_tpl_id,))

    # Insert a new bam_tpl with qualifier 'xenocp'
    query = """insert into bam_tpl (bam_tpl_id, qualifier, status, sample_target_project_id)
            values (nextval('blt_id_seq'), 'xenocp', 'Default', %s) returning bam_tpl_id;"""
    bam_tpl_id_xenocp = raptr.fetch_item_or_fail(query, (stp_id,))

    # Associate loadables with new bam_tpl.
    query = """insert into bam_tpl_read_group values (%s, %s)"""
    for rgid in rg_ids:
        raptr.execute(query, (bam_tpl_id_xenocp, rgid))

    # Get genome for old bam
    query = """select genome_id from bam where bam_id = %s;"""
    genome_id = raptr.fetch_item_or_fail(query, (bam_id,))
    # Add a new bam
    query = """insert into bam (bam_id, bam_tpl_id, status, notes, genome_id)
            values (nextval('blt_id_seq'), %s, 'Normal', NULL, %s) returning bam_id;"""
    bam_id_xenocp = raptr.fetch_item_or_fail(query, (bam_tpl_id_xenocp, genome_id))
    # Update primary_bam_id of bam_tpl with the qualifier
    query = """update bam_tpl set primary_bam_id = %s where bam_tpl_id = %s;"""
    raptr.execute(query, (bam_id_xenocp, bam_tpl_id_xenocp))

    # Get read_group_ids and insert into bam_read_group
    query = """select rg.read_group_id, rg.pu
            from read_group rg inner join sample_target_project_view stpv
            on rg.sample_target_id = stpv.sample_target_id and rg.project_id = stpv.project_id
            where status = 'Normal' and formal_name = %s and target_name = %s
            and project_name = %s and subproject = %s;"""
    res = raptr.execute_fetch(query, (sample, target, project, subproject))
    query = """insert into bam_read_group (bam_id, read_group_id, rgid)
            values (%s, %s, %s);"""
    for (rg_id, pu) in res:
        raptr.execute(query, (bam_id_xenocp, rg_id, pu))

    query = """select control_bam_id from bam_pair where case_bam_id = %s;"""
    control_bam_ids = raptr.fetch_col_or_silent(query, (bam_id,))
    if control_bam_ids:
        # Insert bam_pairs to pair newly created bam_id with control bam_ids with qualifier 'xenocp'
        query = """insert into bam_pair (bam_pair_id, case_bam_id, control_bam_id, qualifier, status, notes)
                values (nextval('blt_id_seq'), %s, %s, 'xenocp', 'Normal', NULL);"""
        for control_bam_id in control_bam_ids:
            raptr.execute(query, (bam_id_xenocp, control_bam_id))


def add_bam_pair_bam_id(raptr, bam_id, **kwargs):
    """ Add a bam pair for a given bam_id
    Args:
        raptr (obj): a Raptr instance
        bam_id (int): bam ID
        kwargs (dict): keyword argument
    Returns:
        None
    """
    query = """select formal_name, target_name, project_name, subproject from sample_target_project_view inner join
            (select sample_target_project_id from bam_and_tpl where bam_id = %s and bam_status = 'Normal' 
            and legacy = false)
            using (sample_target_project_id);"""
    (sample, target, project, subproject) = raptr.fetch_row_or_fail(query, (bam_id,))
    add_bam_pair_stp(raptr, sample, target, project, subproject)


if __name__ == '__main__':
    main()
