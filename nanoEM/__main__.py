#! /usr/bin/env python

from arg_parser import create_parser

def main():
    p = create_parser()
    args = p.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()