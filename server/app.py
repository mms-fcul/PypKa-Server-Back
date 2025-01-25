import datetime
import logging

from flask_cors import CORS

from database import DB_SESSION, BASE, ENGINE

BASE.metadata.create_all(bind=ENGINE)

logging.basicConfig(filename="server.log", level=logging.DEBUG)

root = logging.getLogger("werkzeug")
root.handlers = logging.getLogger().handlers
root.setLevel(logging.INFO)

from flask import Flask
from routes_direct import direct_routes_bp, direct_routes_severelimit_bp
from routes_website import routes_bp

from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

app = Flask(__name__, static_url_path="/static")
app.config["SECRET_KEY"] = str(datetime.datetime.today())
app.config["CORS_HEADERS"] = "Content-Type"

app.register_blueprint(direct_routes_bp)
app.register_blueprint(direct_routes_severelimit_bp)
app.register_blueprint(routes_bp)

limiter = Limiter(
    get_remote_address,
    app=app,
    default_limits=["500 per hour"],
    storage_uri="memory://",
)
limiter.limit("5000/hour")(direct_routes_bp)
limiter.limit("100/hour")(direct_routes_severelimit_bp)

CORS(app, resources={r"/*": {"origins": "*"}})


@app.teardown_appcontext
def shutdown_session(exception=None):
    DB_SESSION.remove()
