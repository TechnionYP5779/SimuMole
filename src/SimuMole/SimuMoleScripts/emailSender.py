import smtplib
import ssl
import os.path as path
import os
import datetime
from email import encoders
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from py3wetransfer import Py3WeTransfer

'''
For import errors run the following commands in the PyCharm terminal:

pip install python-magic-bin==0.4.14
pip install py3wetransfer
'''

api_key = 'pbCr8EGpoO7RvwC7c2eXz4B3QFjgiZhA6WH15eOd'
sender = 'SimuMoleMailbox@gmail.com'
secret = 'Team6!SM'
smtp_server = 'smtp.gmail.com'
port = 465
megabyte = 1000000
size_limit = 10 * megabyte


def send_email(receiver, file_name):
    file_name = "./media/files/" + file_name
    size = os.stat(file_name).st_size
    if size < size_limit:
        send_file_attached(receiver, file_name)
    else:
        send_file_link(receiver, file_name)


def send_file_attached(receiver, file_name):
    message = MIMEMultipart("alternative")
    message["Subject"] = "SimuMole: the files you requested are ready"
    message["From"] = sender
    message["To"] = receiver
    message.attach(MIMEText("The files you requested are attached to this email "
                            "and you can download them at any time.\n\n"
                            "Thank you for using SimuMole.", "plain"))

    attachment = MIMEBase('application', "octet-stream")
    attachment.set_payload(open(file_name, "rb").read())
    encoders.encode_base64(attachment)
    attachment.add_header('Content-Disposition', 'attachment; filename="' + path.basename(file_name) + '"')
    message.attach(attachment)

    context = ssl.create_default_context()

    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender, secret)
        server.sendmail(sender, receiver, message.as_string())


def send_file_link(receiver, file_name):
    x = Py3WeTransfer(api_key)
    url = x.upload_file(file_name, "SimuMole")
    expiration_date = (datetime.datetime.now() + datetime.timedelta(days=7)).strftime('%Y-%m-%d %H:00:00')

    message = MIMEMultipart("alternative")
    message["Subject"] = "SimuMole: the files you requested are ready"
    message["From"] = sender
    message["To"] = receiver
    message.attach(MIMEText("The files you requested are too large to be sent through email.\n"
                            "Thus you can download the files from the following link: " + url + "\n"
                                                                                                "The link is only valid until " + expiration_date +
                            ", so please download your files before then.\n\n"
                            "Thank you for using SimuMole.", "plain"))

    context = ssl.create_default_context()

    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender, secret)
        server.sendmail(sender, receiver, message.as_string())
