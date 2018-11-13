import json
import argparse
import urllib2
import time
import sys
import csv
import getpass
import os
import os.path
import re
import ssl


version = 20160120.1


args = None


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--user', '-u', required=True,
                        help='specify your STARS account')
    parser.add_argument('--release-version', '-r', choices='dr1 dr_early'.split(), required=True,
                        help='specify release version')
    parser.add_argument('--delete-job', '-D', action='store_true',
                        help='delete the job you submitted after your downloading')
    parser.add_argument('--format', '-f', dest='out_format', default='csv', choices=['csv', 'csv.gz', 'sqlite3', 'fits'],
                        help='specify output format')
    parser.add_argument('--nomail', '-M', action='store_true',
                        help='suppress email notice')
    parser.add_argument('--password-env', default='HSC_SSP_CAS_PASSWORD',
                        help='specify the environment variable that has STARS password as its content')
    parser.add_argument('--preview', '-p', action='store_true',
                        help='quick mode (short timeout)')
    parser.add_argument('--skip-syntax-check', '-S', action='store_true',
                        help='skip syntax check')
    parser.add_argument('--api-url', default='https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/',
                        help='for developers')
    parser.add_argument('sql-file', type=argparse.FileType('r'),
                        help='SQL file')

    global args
    args = parser.parse_args()

    credential = {'account_name': args.user, 'password': getPassword()}
    sql = args.__dict__['sql-file'].read()

    job = None

    try:
        if args.preview:
            preview(credential, sql, sys.stdout)
        else:
            job = submitJob(credential, sql, args.out_format)
            blockUntilJobFinishes(credential, job['id'])
            download(credential, job['id'], sys.stdout)
            if args.delete_job:
                deleteJob(credential, job['id'])
    except urllib2.HTTPError, e:
        if e.code == 401:
            print >> sys.stderr, 'invalid id or password.'
        if e.code == 406:
            print >> sys.stderr, e.read()
        else:
            print >> sys.stderr, e
    except QueryError, e:
        print >> sys.stderr, e
    except KeyboardInterrupt:
        if job is not None:
            jobCancel(credential, job['id'])
        raise
    else:
        sys.exit(0)

    sys.exit(1)


class QueryError(Exception):
    pass


def httpJsonPost(url, data):
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})


def httpPost(url, postData, headers):
    req = urllib2.Request(url, postData, headers)
    skipVerifying = None
    try:
        skipVerifying = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
    except AttributeError:
        pass
    if skipVerifying:
        res = urllib2.urlopen(req, context=skipVerifying)
    else:
        res = urllib2.urlopen(req)
    return res


def submitJob(credential, sql, out_format):
    url = args.api_url + 'submit'
    catalog_job = {
        'sql'                     : sql,
        'out_format'              : out_format,
        'include_metainfo_to_body': True,
        'release_version'         : args.release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job, 'nomail': args.nomail, 'skip_syntax_check': args.skip_syntax_check}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def jobStatus(credential, job_id):
    url = args.api_url + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def jobCancel(credential, job_id):
    url = args.api_url + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def preview(credential, sql, out):
    url = args.api_url + 'preview'
    catalog_job = {
        'sql'             : sql,
        'release_version' : args.release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job}
    res = httpJsonPost(url, postData)
    result = json.load(res)

    writer = csv.writer(out)
    # writer.writerow(result['result']['fields'])
    for row in result['result']['rows']:
        writer.writerow(row)

    if result['result']['count'] > len(result['result']['rows']):
        raise QueryError, 'only top %d records are displayed !' % len(result['result']['rows'])


def blockUntilJobFinishes(credential, job_id):
    max_interval = 5 * 60 # sec.
    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError, 'query error: ' + job['error']
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > max_interval:
            interval = max_interval


def download(credential, job_id, out):
    url = args.api_url + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    bufSize = 64 * 1<<10 # 64k
    while True:
        buf = res.read(bufSize)
        out.write(buf)
        if len(buf) < bufSize:
            break


def deleteJob(credential, job_id):
    url = args.api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def getPassword():

    password_from_envvar = os.environ.get(args.password_env, '')
    if password_from_envvar != '':
        return password_from_envvar
    else:
        return getpass.getpass('password? ')


if __name__ == '__main__':
    main()
