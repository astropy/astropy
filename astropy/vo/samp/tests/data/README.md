About
=====

This directory contains test files needed by `astropy.vo.samp`. In particular,
a test SSL certificate and private key are included.

**The SSL certificate and private key should not be used for any other purpose
than for the tests, since the private key is visible to all.**

The key and certificate were generated with:

    $ openssl genrsa -des3 -passout pass:x -out test.pass.key 512
    Generating RSA private key, 512 bit long modulus
    ....++++++++++++
    .....++++++++++++
    e is 65537 (0x10001)

    $ openssl rsa -passin pass:x -in test.pass.key -out test.key
    writing RSA key

    $ openssl req -new -key test.key -out test.csr
    You are about to be asked to enter information that will be incorporated
    into your certificate request.
    What you are about to enter is what is called a Distinguished Name or a DN.
    There are quite a few fields but you can leave some blank
    For some fields there will be a default value,
    If you enter '.', the field will be left blank.
    -----
    Country Name (2 letter code) [AU]:.
    State or Province Name (full name) [Some-State]:.
    Locality Name (eg, city) []:.
    Organization Name (eg, company) [Internet Widgits Pty Ltd]:The Astropy Collaboration
    Organizational Unit Name (eg, section) []:SAMP
    Common Name (e.g. server FQDN or YOUR name) []:
    Email Address []:

    Please enter the following 'extra' attributes
    to be sent with your certificate request
    A challenge password []:
    An optional company name []:

    $ openssl x509 -req -days 365250 -in test.csr -signkey test.key -out test.crt
    Signature ok
    subject=/O=The Astropy Collaboration/OU=SAMP
    Getting Private key

    $ rm test.pass.key test.csr

The certificate expiry date has been set 1000 years in the future, which should
be sufficient for our purposes. The commands above were adapted from the
instructions at http://www.akadia.com/services/ssh_test_certificate.html

The above commands were run twice to produce two certificate/key pairs (since
two were needed for some tests).